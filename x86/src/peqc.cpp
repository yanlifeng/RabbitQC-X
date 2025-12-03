// Created by ylf9811 on 2021/7/6.
//
//

#include "peqc.h"
using namespace std;

#ifdef USE_MPI_IO
#include <mpi.h>
#include <sys/stat.h>

void SkipToLineEndPE(char *data_, int64_t &pos_, const int64_t size_) {
    ASSERT(pos_ < size_);
    while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_) ++pos_;
    if (data_[pos_] == '\r' && pos_ < size_) {
        if (data_[pos_ + 1] == '\n') {
            ++pos_;
        }
    }
}

int64_t GetNextFastqPE(char *data_, int64_t pos_, const int64_t size_) {
    SkipToLineEndPE(data_, pos_, size_);
    ++pos_;

    // find beginning of the next record
    while (data_[pos_] != '@') {
        SkipToLineEndPE(data_, pos_, size_);
        ++pos_;
    }
    int64_t pos0 = pos_;

    SkipToLineEndPE(data_, pos_, size_);
    ++pos_;

    if (data_[pos_] == '@')// previous one was a quality field
        return pos_;
    
    SkipToLineEndPE(data_, pos_, size_);
    ++pos_;
    ASSERT(data_[pos_] == '+');// pos0 was the start of tag
    return pos0;
}
#endif

/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
#ifdef USE_MPI_IO
PeQc::PeQc(CmdInfo *cmd_info1, int my_rank, int comm_size) {
    cmd_info_ = cmd_info1;
    my_rank_ = my_rank;
    comm_size_ = comm_size;
    
    // Initialize basic variables
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    nowChunkId = 1;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != string::npos;

    // Check for unsupported MPI configurations
    if (in_is_zip_) {
        fprintf(stderr, "Error: MPI distributed reading does not support compressed (.gz) input files.\n");
        fprintf(stderr, "Please use uncompressed FASTQ files for MPI mode.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (!cmd_info1->splitWrite_ && out_is_zip_) {
        fprintf(stderr, "Error: MPI single file write mode does not support compressed (.gz) output.\n");
        fprintf(stderr, "Please use --splitWrite for compressed output, or use uncompressed output for single file mode.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Check if both input files have the same size
    struct stat st1, st2;
    if (stat(cmd_info1->in_file_name1_.c_str(), &st1) != 0) {
        fprintf(stderr, "Error: cannot stat file %s\n", cmd_info1->in_file_name1_.c_str());
        exit(1);
    }
    if (stat(cmd_info1->in_file_name2_.c_str(), &st2) != 0) {
        fprintf(stderr, "Error: cannot stat file %s\n", cmd_info1->in_file_name2_.c_str());
        exit(1);
    }
    
    int64_t real_file_size1 = st1.st_size;
    int64_t real_file_size2 = st2.st_size;
    
    if (!in_is_zip_ && real_file_size1 != real_file_size2) {
        fprintf(stderr, "Error: PE input files have different sizes: %lld vs %lld\n", 
                real_file_size1, real_file_size2);
        exit(1);
    }

    // Calculate initial approximate division (based on file1)
    int64_t pre_size = (real_file_size1 + comm_size - 1) / comm_size;
    int64_t start_pos = pre_size * my_rank;
    int64_t end_pos = start_pos + pre_size;
    if (end_pos > real_file_size1) end_pos = real_file_size1;

    // Adjust boundaries to FASTQ record boundaries (only check file1)
    int64_t right_pos = 0;
    if (!in_is_zip_) {
        FILE *pre_fp = fopen(cmd_info1->in_file_name1_.c_str(), "rb");
        if (pre_fp == NULL) {
            fprintf(stderr, "Error: cannot open file %s\n", cmd_info1->in_file_name1_.c_str());
            exit(1);
        }
        fseek(pre_fp, start_pos, SEEK_SET);
        char *tmp_chunk = new char[1 << 20];
        int res_size = fread(tmp_chunk, sizeof(char), 1 << 20, pre_fp);
        
        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0) {
            right_pos = 0;
        } else {
            right_pos = GetNextFastqPE(tmp_chunk, 0, res_size);
        }
        fclose(pre_fp);
        delete[] tmp_chunk;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    int64_t now_pos = right_pos + start_pos;
    int64_t *now_poss = new int64_t[comm_size];
    now_poss[0] = now_pos;
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank) {
        MPI_Send(&now_pos, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    } else {
        for (int ii = 1; ii < comm_size; ii++) {
            MPI_Recv(&(now_poss[ii]), 1, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) {
        for (int ii = 1; ii < comm_size; ii++) {
            MPI_Send(now_poss, comm_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(now_poss, comm_size, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start_pos_ = now_poss[my_rank];
    if (my_rank == comm_size - 1) {
        end_pos_ = real_file_size1;
    } else {
        end_pos_ = now_poss[my_rank + 1];
    }
    
    delete[] now_poss;

    debug_printf("MPI PE file splitting: rank %d, start: %lld, end: %lld\n", 
                 my_rank, start_pos_, end_pos_);

    // Initialize MPI RMA windows for single file write mode
    mpi_win1_ = MPI_WIN_NULL;
    mpi_win2_ = MPI_WIN_NULL;
    global_offset_ptr1_ = NULL;
    global_offset_ptr2_ = NULL;
    mpi_file1_ = MPI_FILE_NULL;
    mpi_file2_ = MPI_FILE_NULL;
    
    if (!cmd_info1->splitWrite_ && cmd_info1->write_data_ && !out_is_zip_) {
        // Allocate shared memory for global offset (only on rank 0)
        MPI_Aint size = (my_rank == 0) ? sizeof(long long) : 0;
        int disp_unit = sizeof(long long);
        
        // For file1
        MPI_Win_allocate(size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, 
                         &global_offset_ptr1_, &mpi_win1_);
        if (my_rank == 0) {
            *global_offset_ptr1_ = 0;
        }
        MPI_Win_lock_all(0, mpi_win1_);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name1_.c_str(),
                      MPI_MODE_CREATE | MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &mpi_file1_);
        
        // For file2 (if not interleaved output)
        if (cmd_info1->interleaved_out_ == 0 && cmd_info1->out_file_name2_.length() > 0) {
            MPI_Win_allocate(size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, 
                             &global_offset_ptr2_, &mpi_win2_);
            if (my_rank == 0) {
                *global_offset_ptr2_ = 0;
            }
            MPI_Win_lock_all(0, mpi_win2_);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name2_.c_str(),
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &mpi_file2_);
        }
        
        debug_printf("Rank %d: MPI PE single file write initialized\n", my_rank);
    }

    // Continue with common initialization
    out_queue1_ = NULL;
    out_queue2_ = NULL;
#else
PeQc::PeQc(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    out_queue1_ = NULL;
    out_queue2_ = NULL;
    nowChunkId = 1;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != string::npos;
#endif
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    out_queue1_ = NULL;
    out_queue2_ = NULL;
    nowChunkId = 1;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != string::npos;
    if (cmd_info1->write_data_) {
        out_queue1_ = new CIPair[1 << 20];
        queue1P1 = 0;
        queue1P2 = 0;
        queueNumNow1 = 0;
        queueSizeLim1 = 1 << 5;
        if (cmd_info_->interleaved_out_ == 0) {
            out_queue2_ = new CIPair[1 << 20];
            queue2P1 = 0;
            queue2P2 = 0;
            queueNumNow2 = 0;
            queueSizeLim2 = 1 << 5;
        }
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                pigzQueueNumNow1 = 0;
                pigzQueueSizeLim1 = 1 << 5;
                string out_name1 = cmd_info1->out_file_name1_;
                out_name1 = out_name1.substr(0, out_name1.find(".gz"));
                out_stream1_.open(out_name1);
                out_stream1_.close();

                pigzQueueNumNow2 = 0;
                pigzQueueSizeLim2 = 1 << 5;
                string out_name2 = cmd_info1->out_file_name2_;
                out_name2 = out_name2.substr(0, out_name2.find(".gz"));
                out_stream2_.open(out_name2);
                out_stream2_.close();
                debug_printf("now use pigz to compress output data\n");

            } else {
                debug_printf("open gzip stream1 %s\n", cmd_info1->out_file_name1_.c_str());
                debug_printf("open gzip stream2 %s\n", cmd_info1->out_file_name2_.c_str());
                zip_out_stream1 = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
                gzsetparams(zip_out_stream1, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream1, 1024 * 1024);
                zip_out_stream2 = gzopen(cmd_info1->out_file_name2_.c_str(), "w");
                gzsetparams(zip_out_stream2, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream2, 1024 * 1024);
            }
        } else {
#ifdef USE_MPI_IO
            // Don't open regular file stream if using MPI single file mode
            if (cmd_info1->splitWrite_) {
#endif
                debug_printf("open stream1 %s\n", cmd_info1->out_file_name1_.c_str());
                if (cmd_info_->interleaved_out_ == 0)
                    debug_printf("open stream2 %s\n", cmd_info1->out_file_name2_.c_str());
                out_stream1_.open(cmd_info1->out_file_name1_);
                if (cmd_info_->interleaved_out_ == 0){
                    out_stream2_.open(cmd_info1->out_file_name2_);
                }
#ifdef USE_MPI_IO
            }
#endif
        }
    }
    duplicate_ = NULL;
    if (cmd_info1->state_duplicate_) {
        duplicate_ = new Duplicate(cmd_info1);
    }
    umier_ = NULL;
    if (cmd_info1->add_umi_) {
        umier_ = new Umier(cmd_info1);
    }
    if (cmd_info1->use_pugz_) {
        pugzQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 5);
        pugzQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 5);
    }
    if (cmd_info1->use_pigz_) {
        pigzQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>;
        pigzLast1.first = new char[1 << 24];
        pigzLast1.second = 0;
        pigzQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>;
        pigzLast2.first = new char[1 << 24];
        pigzLast2.second = 0;
    }
    if(cmd_info1->do_correction_with_care_) {
        careQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 20);
        careQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 20);
    }
    changeNum = 0;
    careStartWrite = 0;
    careDone1 = 0;
    careDone2 = 0;
    pugzDone1 = 0;
    pugzDone2 = 0;
    producerDone = 0;
    writerDone1 = 0;
    writerDone2 = 0;
}


PeQc::~PeQc() {
    delete filter_;
    if (cmd_info_->write_data_) {
        delete out_queue1_;
        if (cmd_info_->interleaved_out_ == 0)
            delete out_queue2_;
    }
    if (cmd_info_->state_duplicate_) {
        delete duplicate_;
    }
    if (cmd_info_->add_umi_) {
        delete umier_;
    }
}

void PeQc::careProcess() {
    
    //fprintf(stderr, "pairmode %s\n", cmd_info_->pairmode_.c_str());
    //fprintf(stderr, "coverage %d\n", cmd_info_->coverage_);
    vector<char*>paras;
    paras.push_back("./RabbitQCPlus");
    paras.push_back("-i");
    paras.push_back((char*)(cmd_info_->in_file_name1_.data()));
    paras.push_back("-i");
    paras.push_back((char*)(cmd_info_->in_file_name2_.data()));

    paras.push_back("-d");
    paras.push_back("./");

    paras.push_back("-o");
    paras.push_back("tmp1.fq");
    
    paras.push_back("-o");
    paras.push_back("tmp2.fq");

    paras.push_back("-c");
    string str_coverage = to_string(cmd_info_->coverage_);
    //fprintf(stderr, "str_coverage %s\n", str_coverage.c_str());
    paras.push_back((char*)(str_coverage.data()));

    paras.push_back("-t");
    string str_thread = to_string(cmd_info_->correct_threadnum_);
    //fprintf(stderr, "str_thread %s\n", str_thread.c_str());
    paras.push_back((char*)(str_thread.data()));

    paras.push_back("--pairmode");
    if(cmd_info_->pairmode_ == "SE") paras.push_back("SE");
    else paras.push_back("PE");


    string str_hashmaps = to_string(cmd_info_->hashmaps_);
    if(cmd_info_->hashmaps_ != 48) {
        paras.push_back("--hashmaps");
        paras.push_back((char*)(str_hashmaps.data()));
    }
    string str_kmerlength = to_string(cmd_info_->kmerlength_);
    if(cmd_info_->kmerlength_ != 0) {
        paras.push_back("--kmerlength");
        paras.push_back((char*)(str_kmerlength.data()));
    }
    if(cmd_info_->enforceHashmapCount_ != false) {
        paras.push_back("--enforceHashmapCount");
    }
    if(cmd_info_->useQualityScores_ != false) {
        paras.push_back("--useQualityScores");
    }
    string str_qualityScoreBits = to_string(cmd_info_->qualityScoreBits_);
    if(cmd_info_->qualityScoreBits_ != 8) {
        paras.push_back("--qualityScoreBits");
        paras.push_back((char*)(str_qualityScoreBits.data()));
    }
    if(cmd_info_->excludeAmbiguous_ != false) {
        paras.push_back("--excludeAmbiguous");
    }
    string str_maxmismatchratio = to_string(cmd_info_->maxmismatchratio_);
    if(cmd_info_->maxmismatchratio_ != 0.200000) {
        paras.push_back("--maxmismatchratio");
        paras.push_back((char*)(str_maxmismatchratio.data()));
    }
    string str_minalignmentoverlap = to_string(cmd_info_->minalignmentoverlap_);
    if(cmd_info_->minalignmentoverlap_ != 30) {
        paras.push_back("--minalignmentoverlap");
        paras.push_back((char*)(str_minalignmentoverlap.data()));
    }
    string str_minalignmentoverlapratio = to_string(cmd_info_->minalignmentoverlapratio_);
    if(cmd_info_->minalignmentoverlapratio_ != 0.300000) {
        paras.push_back("--minalignmentoverlapratio");
        paras.push_back((char*)(str_minalignmentoverlapratio.data()));
    }
    string str_errorfactortuning = to_string(cmd_info_->errorfactortuning_);
    if(cmd_info_->errorfactortuning_ != 0.060000) {
        paras.push_back("--errorfactortuning");
        paras.push_back((char*)(str_errorfactortuning.data()));
    }
    string str_coveragefactortuning = to_string(cmd_info_->coveragefactortuning_);
    if(cmd_info_->coveragefactortuning_ != 0.600000) {
        paras.push_back("--coveragefactortuning");
        paras.push_back((char*)(str_coveragefactortuning.data()));
    }
    if(cmd_info_->showProgress_ != false) {
        paras.push_back("--showProgress");
    }
    string str_tempdir = cmd_info_->tempdir_;
    if(cmd_info_->tempdir_ != "") {
        paras.push_back("--tempdir");
        paras.push_back((char*)(str_tempdir.data()));
    }
    string str_save_preprocessedreads_to = cmd_info_->save_preprocessedreads_to_;
    if(cmd_info_->save_preprocessedreads_to_ != "") {
        paras.push_back("--save_preprocessedreads_to");
        paras.push_back((char*)(str_save_preprocessedreads_to.data()));
    }
    string str_load_preprocessedreads_from = cmd_info_->load_preprocessedreads_from_;
    if(cmd_info_->load_preprocessedreads_from_ != "") {
        paras.push_back("--load_preprocessedreads_from");
        paras.push_back((char*)(str_load_preprocessedreads_from.data()));
    }
    string str_save_hashtables_to = cmd_info_->save_hashtables_to_;
    if(cmd_info_->save_hashtables_to_ != "") {
        paras.push_back("--save_hashtables_to");
        paras.push_back((char*)(str_save_hashtables_to.data()));
    }
    string str_load_hashtables_from = cmd_info_->load_hashtables_from_;
    if(cmd_info_->load_hashtables_from_ != "") {
        paras.push_back("--load_hashtables_from");
        paras.push_back((char*)(str_load_hashtables_from.data()));
    }
    string str_memHashtables = cmd_info_->memHashtables_;
    if(cmd_info_->memHashtables_ != "") {
        paras.push_back("--memHashtables");
        paras.push_back((char*)(str_memHashtables.data()));
    }
    string str_memTotal = cmd_info_->memTotal_;
    if(cmd_info_->memTotal_ != "") {
        paras.push_back("--memTotal");
        paras.push_back((char*)(str_memTotal.data()));
    }
    string str_hashloadfactor = to_string(cmd_info_->hashloadfactor_);
    if(cmd_info_->hashloadfactor_ != 0.800000) {
        paras.push_back("--hashloadfactor");
        paras.push_back((char*)(str_hashloadfactor.data()));
    }
    string str_fixedNumberOfReads = to_string(cmd_info_->fixedNumberOfReads_);
    if(cmd_info_->fixedNumberOfReads_ != 0) {
        paras.push_back("--fixedNumberOfReads");
        paras.push_back((char*)(str_fixedNumberOfReads.data()));
    }
    if(cmd_info_-> singlehash_!= false) {
        paras.push_back("--singlehash");
    }
    if(cmd_info_->correctionQualityLabels_ != false) {
        paras.push_back("--correctionQualityLabels");
    }
    if(cmd_info_->candidateCorrection_ != false) {
        paras.push_back("--candidateCorrection");
    }
    string str_candidateCorrectionNewColumns = to_string(cmd_info_->candidateCorrectionNewColumns_);
    if(cmd_info_->candidateCorrectionNewColumns_ != 15) {
        paras.push_back("--candidateCorrectionNewColumns");
        paras.push_back((char*)(str_candidateCorrectionNewColumns.data()));
    }
    string str_correctionType = to_string(cmd_info_->correctionType_);
    if(cmd_info_->correctionType_ != 0) {
        paras.push_back("--correctionType");
        paras.push_back((char*)(str_correctionType.data()));
    }
    string str_correctionTypeCands = to_string(cmd_info_->correctionTypeCands_);
    if(cmd_info_->correctionTypeCands_ != 0) {
        paras.push_back("--correctionTypeCands");
        paras.push_back((char*)(str_correctionTypeCands.data()));
    }






    //for(int i = 0; i < paras.size(); i++) {
    //    fprintf(stderr, "%s ", paras[i]);
    //}
    //fprintf(stderr, "\n");

    fprintf(stderr, "start care part...\n");

    //fprintf(stderr, "now output to queue, %p %p %p\n", careQueue1, careQueue2, &producerDone);
    main_correction(paras.size(), &(paras[0]), careQueue1, careQueue2, &producerDone, &careStartWrite, &changeNum);

    fprintf(stderr, "care end\n");
    //fprintf(stderr, "care1 queue size %d\n", careQueue1->size_approx());
    //fprintf(stderr, "care2 queue size %d\n", careQueue2->size_approx());
    //fprintf(stderr, "care change size %d\n", changeNum);

    careDone1 = 1;
    careDone2 = 1;

}



string PeQc::Read2String(neoReference &ref) {
    return string((char *) ref.base + ref.pname, ref.lname) + "\n" +
        string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
        string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
        string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
}

void PeQc::Read2Chars(neoReference &ref, char *out_data, int &pos) {
    memcpy(out_data + pos, ref.base + ref.pname, ref.lname);
    pos += ref.lname;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pseq, ref.lseq);
    pos += ref.lseq;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pstrand, ref.lstrand);
    pos += ref.lstrand;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pqual, ref.lqual);
    pos += ref.lqual;
    out_data[pos++] = '\n';
}

void PeQc::ProducerPeInterFastqTask(string file, rabbit::fq::FastqDataPool *fastq_data_pool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    debug_printf("ProducerPeInterFastqTask started\n");
    double t0 = GetTime();
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastq_data_pool, "", in_is_zip_);
    int64_t n_chunks = 0;
    while (true) {
        rabbit::fq::FastqDataChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextInterChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        dq.Push(n_chunks, fqdatachunk);
    }

    dq.SetCompleted();
    delete fqFileReader;
    //cout << "file " << file << " has " << n_chunks << " chunks" << endl;
    debug_printf("producer cost %.3f\n", GetTime() - t0);
}


void PeQc::ProducerPeFastqTask(string file, string file2, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    debug_printf("ProducerPeFastqTask started\n");
    double t0 = GetTime();
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
#ifdef USE_MPI_IO
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, file2, in_is_zip_, tmpSize, start_pos_, end_pos_);
#else
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, file2, in_is_zip_, tmpSize);
#endif
    int n_chunks = 0;


    if (cmd_info_->use_pugz_) {
        pair<char *, int> last1;
        pair<char *, int> last2;
        last1.first = new char[1 << 20];
        last1.second = 0;
        last2.first = new char[1 << 20];
        last2.second = 0;
        while (true) {
            rabbit::fq::FastqDataPairChunk *fqdatachunk;
            //fqdatachunk = fqFileReader->readNextPairChunkParallel(pugzQueue1, pugzQueue2, &pugzDone1, &pugzDone2, last1,
            fqdatachunk = fqFileReader->readNextPairChunk(pugzQueue1, pugzQueue2, &pugzDone1, &pugzDone2, last1,
                    last2);
            if (fqdatachunk == NULL) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
        delete[] last1.first;
        delete[] last2.first;
    } else if(cmd_info_->do_correction_with_care_){
        pair<char *, int> last1;
        pair<char *, int> last2;
        last1.first = new char[1 << 20];
        last1.second = 0;
        last2.first = new char[1 << 20];
        last2.second = 0;
        while (true) {
            rabbit::fq::FastqDataPairChunk *fqdatachunk;
            //fqdatachunk = fqFileReader->readNextPairChunkParallel(careQueue1, careQueue2, &careDone1, &careDone2, last1,
            fqdatachunk = fqFileReader->readNextPairChunk(careQueue1, careQueue2, &careDone1, &careDone2, last1,
                    last2);
            if (fqdatachunk == NULL) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
        delete[] last1.first;
        delete[] last2.first;
    } else {
        long long tot_read_size = 0;
        long long tot_read_size2 = 0;
        while (true) {
            rabbit::fq::FastqDataPairChunk *fqdatachunk;
            //fqdatachunk = fqFileReader->readNextPairChunkParallel();
            fqdatachunk = fqFileReader->readNextPairChunk();
            if (fqdatachunk == NULL) break;
            tot_read_size += fqdatachunk->left_part->size;
            tot_read_size2 += fqdatachunk->right_part->size;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
        // printf("tot_read_size: %lld\n", tot_read_size);
        // printf("tot_read_size2: %lld\n", tot_read_size2);
        // printf("n_chunks: %lld\n", n_chunks);
    }


    dq.SetCompleted();
    delete fqFileReader;
    //cout << "file " << file << " has " << n_chunks << " chunks" << endl;
    debug_printf("producer cost %.3f\n", GetTime() - t0);
}

/**
 * @brief get fastq data chunks from the data queue and do QC for them
 * @param thread_info : thread information
 * @param fastq_data_pool :a fastq data pool, it will be used to release data chunk
 * @param dq : data queue
 */
void PeQc::ConsumerPeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
    while (dq.Pop(id, fqdatachunk)) {
        vector<neoReference> data1, data2;
        vector<neoReference> pass_data1, pass_data2;
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->left_part), data1, true);
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->right_part), data2, true);
        ASSERT(data1.size() == data2.size());
        int out_len1 = 0, out_len2 = 0;
        int b_size = min(data1.size(), data2.size());
        for (int i = 0; i < b_size; i++) {
            auto item1 = data1[i];
            auto item2 = data2[i];
            thread_info->pre_state1_->StateInfo(item1);
            thread_info->pre_state2_->StateInfo(item2);
            if (cmd_info_->state_duplicate_) {
                duplicate_->statPair(item1, item2);
            }
            if (cmd_info_->add_umi_) {
                umier_->ProcessPe(item1, item2);
            }

            //do pe sequence trim
            bool trim_res1 = filter_->TrimSeq(item1, cmd_info_->trim_front1_, cmd_info_->trim_tail1_);
            bool trim_res2 = filter_->TrimSeq(item2, cmd_info_->trim_front2_, cmd_info_->trim_tail2_);

            if (trim_res1 && trim_res2 && cmd_info_->trim_polyg_) {
                PolyX::trimPolyG(item1, item2, cmd_info_->trim_poly_len_);
            }

            if (trim_res1 && trim_res2 && cmd_info_->trim_polyx_) {
                PolyX::trimPolyX(item1, item2, cmd_info_->trim_poly_len_);
            }
            //do pe overlap analyze
            OverlapRes overlap_res;
            if (trim_res1 && trim_res2 && cmd_info_->analyze_overlap_) {
                overlap_res = Adapter::AnalyzeOverlap(item1, item2, cmd_info_->overlap_diff_limit_,
                        cmd_info_->overlap_require_);
                int now_size = cmd_info_->max_insert_size_;
                if (overlap_res.overlapped) {
                    if (overlap_res.offset > 0)
                        now_size = item1.lseq + item2.lseq - overlap_res.overlap_len;
                    else
                        now_size = overlap_res.overlap_len;
                }
                now_size = min(now_size, cmd_info_->max_insert_size_);
                thread_info->insert_size_dist_[now_size]++;
            }
            if (trim_res1 && trim_res2 && cmd_info_->correct_data_) {
                Adapter::CorrectData(item1, item2, overlap_res, cmd_info_->isPhred64_);
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
                int trimmed= false;
                if (cmd_info_->print_what_trimmed_) {
                    trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len,
                            thread_info->aft_state1_->adapter_map_,
                            thread_info->aft_state2_->adapter_map_, cmd_info_->adapter_len_lim_);
                } else {
                    trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len);
                }
                if (trimmed) {
                    thread_info->aft_state1_->AddTrimAdapter();
                    thread_info->aft_state1_->AddTrimAdapter();
                    thread_info->aft_state1_->AddTrimAdapterBase(trimmed);
                }
                int res1 = 0, res2 = 0;
                bool is_trimmed1 = trimmed;
                bool is_trimmed2 = trimmed;

                if (!trimmed) {
                    if (cmd_info_->detect_adapter1_) {
                        int res1;
                        if (cmd_info_->print_what_trimmed_) {
                            res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_,
                                    thread_info->aft_state1_->adapter_map_,
                                    cmd_info_->adapter_len_lim_, false);
                        } else {
                            res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_, false);
                        }
                        if (res1) {
                            is_trimmed1 = true;
                            thread_info->aft_state1_->AddTrimAdapter();
                            thread_info->aft_state1_->AddTrimAdapterBase(res1);
                        }
                    }
                    if (cmd_info_->detect_adapter2_) {
                        int res2;
                        if (cmd_info_->print_what_trimmed_) {
                            res2 = Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_,
                                    thread_info->aft_state2_->adapter_map_,
                                    cmd_info_->adapter_len_lim_, true);
                        } else {
                            res2 = Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_, true);
                        }
                        if (res2) {
                            is_trimmed2 = true;
                            thread_info->aft_state1_->AddTrimAdapter();
                            thread_info->aft_state1_->AddTrimAdapterBase(res2);
                        }

                    }


                }

                if(cmd_info_->adapter_from_fasta_.size() > 0) {
                    if (cmd_info_->print_what_trimmed_) {
                        res1 = Adapter::TrimAdapters(item1, cmd_info_->adapter_from_fasta_,
                                thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                    } else {
                        res1 = Adapter::TrimAdapters(item1, cmd_info_->adapter_from_fasta_, false);
                    }
                    if (res1) {
                        if(!is_trimmed1) thread_info->aft_state1_->AddTrimAdapter();
                        thread_info->aft_state1_->AddTrimAdapterBase(res1);
                    }

                    if (cmd_info_->print_what_trimmed_) {
                        res2 = Adapter::TrimAdapters(item2, cmd_info_->adapter_from_fasta_,
                                thread_info->aft_state2_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                    } else {
                        res2 = Adapter::TrimAdapters(item2, cmd_info_->adapter_from_fasta_, false);
                    }
                    if (res2) {
                        if(!is_trimmed2) thread_info->aft_state1_->AddTrimAdapter();
                        thread_info->aft_state1_->AddTrimAdapterBase(res2);
                    }
                }

            }


            //do filer in refs
            int filter_res1 = filter_->ReadFiltering(item1, trim_res1, cmd_info_->isPhred64_);
            int filter_res2 = filter_->ReadFiltering(item2, trim_res2, cmd_info_->isPhred64_);

            int filter_res = max(filter_res1, filter_res2);


            if (filter_res == 0) {
                thread_info->aft_state1_->StateInfo(item1);
                thread_info->aft_state2_->StateInfo(item2);
                thread_info->aft_state1_->AddPassReads();
                thread_info->aft_state1_->AddPassReads();
                if (cmd_info_->write_data_) {
                    pass_data1.push_back(item1);
                    pass_data2.push_back(item2);
                    out_len1 += item1.lname + item1.lseq + item1.lstrand + item1.lqual + 4;
                    out_len2 += item2.lname + item2.lseq + item2.lstrand + item2.lqual + 4;
                }
            } else if (filter_res == 2) {
                thread_info->aft_state1_->AddFailShort();
                thread_info->aft_state1_->AddFailShort();
            } else if (filter_res == 3) {
                thread_info->aft_state1_->AddFailLong();
                thread_info->aft_state1_->AddFailLong();
            } else if (filter_res == 4) {
                thread_info->aft_state1_->AddFailLowq();
                thread_info->aft_state1_->AddFailLowq();
            } else if (filter_res == 1) {
                thread_info->aft_state1_->AddFailN();
                thread_info->aft_state1_->AddFailN();
            }
        }
        if (cmd_info_->write_data_) {
            if (cmd_info_->interleaved_out_) {
                char *out_data = new char[out_len1 + out_len2];
                int pos = 0;
                int len = min(pass_data1.size(), pass_data2.size());
                for (int i = 0; i < len; i++) {
                    auto item1 = pass_data1[i];
                    auto item2 = pass_data2[i];
                    Read2Chars(item1, out_data, pos);
                    Read2Chars(item2, out_data, pos);
                }

                if(cmd_info_->notKeepOrder_ == 0) {
                    while(nowChunkId != id) {
                        usleep(100);
                    }
                }

                mylock.lock();
                    while (queueNumNow1 >= queueSizeLim1) {
                    //debug_printf("waiting to push a chunk to out queue1\n");
                    usleep(100);
                }
                out_queue1_[queue1P2++] = {out_data, pos};
                queueNumNow1++;
                nowChunkId++;
                mylock.unlock();
            } else {
                if (pass_data1.size() > 0 && pass_data2.size() > 0) {
                    char *out_data1 = new char[out_len1];
                    int pos = 0;
                    for (auto item: pass_data1) {
                        Read2Chars(item, out_data1, pos);
                    }
                    ASSERT(pos == out_len1);
                    char *out_data2 = new char[out_len2];
                    pos = 0;
                    for (auto item: pass_data2) {
                        Read2Chars(item, out_data2, pos);
                    }
                    ASSERT(pos == out_len2);

                    if(cmd_info_->notKeepOrder_ == 0) {
                        while(nowChunkId != id) {
                            usleep(100);
                        }
                    }
                    mylock.lock();
                    while (queueNumNow1 >= queueSizeLim1 || queueNumNow2 >= queueSizeLim2) {
                        //debug_printf("waiting to push a chunk to out queue1\n");
                        usleep(100);
                    }
                    out_queue1_[queue1P2++] = {out_data1, out_len1};
                    queueNumNow1++;
                    out_queue2_[queue2P2++] = {out_data2, out_len2};
                    queueNumNow2++;
                    nowChunkId++;
                    mylock.unlock();
                }
            }
        }
        fastqPool->Release(fqdatachunk->left_part);
        fastqPool->Release(fqdatachunk->right_part);
    }
    done_thread_number_++;
}


void PeQc::ConsumerPeInterFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;
    while (dq.Pop(id, fqdatachunk)) {
        vector<neoReference> data;
        vector<neoReference> pass_data1, pass_data2;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);
        ASSERT(data.size() % 2 == 0);
        int out_len1 = 0, out_len2 = 0;
        for (int i = 0; i + 2 <= data.size(); i += 2) {
            auto item1 = data[i];
            auto item2 = data[i + 1];
            thread_info->pre_state1_->StateInfo(item1);
            thread_info->pre_state2_->StateInfo(item2);
            if (cmd_info_->state_duplicate_) {
                duplicate_->statPair(item1, item2);
            }
            if (cmd_info_->add_umi_) {
                umier_->ProcessPe(item1, item2);
            }

            //do pe sequence trim
            bool trim_res1 = filter_->TrimSeq(item1, cmd_info_->trim_front1_, cmd_info_->trim_tail1_);
            bool trim_res2 = filter_->TrimSeq(item2, cmd_info_->trim_front2_, cmd_info_->trim_tail2_);

            if (trim_res1 && trim_res2 && cmd_info_->trim_polyg_) {
                PolyX::trimPolyG(item1, item2, cmd_info_->trim_poly_len_);
            }

            if (trim_res1 && trim_res2 && cmd_info_->trim_polyx_) {
                PolyX::trimPolyX(item1, item2, cmd_info_->trim_poly_len_);
            }

            //do pe overlap analyze
            OverlapRes overlap_res;
            if (trim_res1 && trim_res2 && cmd_info_->analyze_overlap_) {
                overlap_res = Adapter::AnalyzeOverlap(item1, item2, cmd_info_->overlap_diff_limit_,
                        cmd_info_->overlap_require_);
                int now_size = cmd_info_->max_insert_size_;
                if (overlap_res.overlapped) {
                    if (overlap_res.offset > 0)
                        now_size = item1.lseq + item2.lseq - overlap_res.overlap_len;
                    else
                        now_size = overlap_res.overlap_len;
                }
                now_size = min(now_size, cmd_info_->max_insert_size_);
                thread_info->insert_size_dist_[now_size]++;
            }
            if (trim_res1 && trim_res2 && cmd_info_->correct_data_) {
                Adapter::CorrectData(item1, item2, overlap_res, cmd_info_->isPhred64_);
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
                int trimmed;
                if (cmd_info_->print_what_trimmed_) {
                    trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len,
                            thread_info->aft_state1_->adapter_map_,
                            thread_info->aft_state2_->adapter_map_, cmd_info_->adapter_len_lim_);
                } else {
                    trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len);
                }
                if (trimmed) {
                    thread_info->aft_state1_->AddTrimAdapter();
                    thread_info->aft_state1_->AddTrimAdapter();
                    thread_info->aft_state1_->AddTrimAdapterBase(trimmed);
                }
                if (!trimmed) {
                    int res1, res2;
                    if (cmd_info_->detect_adapter1_) {
                        int res1;
                        if (cmd_info_->print_what_trimmed_) {
                            res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_,
                                    thread_info->aft_state1_->adapter_map_,
                                    cmd_info_->adapter_len_lim_, false);
                        } else {
                            res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_, false);
                        }
                    }
                    if (cmd_info_->detect_adapter2_) {
                        int res2;
                        if (cmd_info_->print_what_trimmed_) {
                            res2 = Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_,
                                    thread_info->aft_state2_->adapter_map_,
                                    cmd_info_->adapter_len_lim_, true);
                        } else {
                            res2 = Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_, true);
                        }
                    }
                    if (res1) {
                        thread_info->aft_state1_->AddTrimAdapter();
                        thread_info->aft_state1_->AddTrimAdapterBase(res1);
                    }
                    if (res2) {
                        thread_info->aft_state1_->AddTrimAdapter();
                        thread_info->aft_state1_->AddTrimAdapterBase(res2);
                    }
                }
            }


            //do filer in refs
            int filter_res1 = filter_->ReadFiltering(item1, trim_res1, cmd_info_->isPhred64_);
            int filter_res2 = filter_->ReadFiltering(item2, trim_res2, cmd_info_->isPhred64_);


            int filter_res = max(filter_res1, filter_res2);


            if (filter_res == 0) {
                thread_info->aft_state1_->StateInfo(item1);
                thread_info->aft_state2_->StateInfo(item2);
                thread_info->aft_state1_->AddPassReads();
                thread_info->aft_state1_->AddPassReads();
                if (cmd_info_->write_data_) {
                    pass_data1.push_back(item1);
                    pass_data2.push_back(item2);
                    out_len1 += item1.lname + item1.lseq + item1.lstrand + item1.lqual + 4;
                    out_len2 += item2.lname + item2.lseq + item2.lstrand + item2.lqual + 4;
                }
            } else if (filter_res == 2) {
                thread_info->aft_state1_->AddFailShort();
                thread_info->aft_state1_->AddFailShort();
            } else if (filter_res == 3) {
                thread_info->aft_state1_->AddFailLong();
                thread_info->aft_state1_->AddFailLong();
            } else if (filter_res == 4) {
                thread_info->aft_state1_->AddFailLowq();
                thread_info->aft_state1_->AddFailLowq();
            } else if (filter_res == 1) {
                thread_info->aft_state1_->AddFailN();
                thread_info->aft_state1_->AddFailN();
            }
        }
        if (cmd_info_->write_data_) {
            if (cmd_info_->interleaved_out_) {
                char *out_data = new char[out_len1 + out_len2];
                int pos = 0;
                int len = min(pass_data1.size(), pass_data2.size());
                for (int i = 0; i < len; i++) {
                    auto item1 = pass_data1[i];
                    auto item2 = pass_data2[i];
                    Read2Chars(item1, out_data, pos);
                    Read2Chars(item2, out_data, pos);
                }
                mylock.lock();
                    while (queueNumNow1 >= queueSizeLim1) {
                    //debug_printf("waiting to push a chunk to out queue1\n");
                    usleep(100);
                }
                out_queue1_[queue1P2++] = {out_data, pos};
                queueNumNow1++;
                mylock.unlock();
            } else {
                if (pass_data1.size() > 0 && pass_data2.size() > 0) {
                    char *out_data1 = new char[out_len1];
                    int pos = 0;
                    for (auto item: pass_data1) {
                        Read2Chars(item, out_data1, pos);
                    }
                    ASSERT(pos == out_len1);
                    char *out_data2 = new char[out_len2];
                    pos = 0;
                    for (auto item: pass_data2) {
                        Read2Chars(item, out_data2, pos);
                    }
                    ASSERT(pos == out_len2);

                    mylock.lock();
                    while (queueNumNow1 >= queueSizeLim1 || queueNumNow2 >= queueSizeLim2) {
                        //debug_printf("waiting to push a chunk to out queue1\n");
                        usleep(100);
                    }
                    out_queue1_[queue1P2++] = {out_data1, out_len1};
                    queueNumNow1++;
                    out_queue2_[queue2P2++] = {out_data2, out_len2};
                    queueNumNow2++;
                    mylock.unlock();
                }
            }
        }

        fastqPool->Release(fqdatachunk);
    }
    done_thread_number_++;
}


/**
 * @brief a function to write pe data from out_data1 queue to file1
 */
void PeQc::WriteSeFastqTask1() {
    debug_printf("WriteSeFastqTask1 started\n");
    double t0 = GetTime();
    int cnt = 0;
    bool overWhile = 0;
    pair<char *, int> now;
    while (true) {
        while (queueNumNow1 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                break;
            }
            usleep(100);
        }
        if (overWhile) break;
        now = out_queue1_[queue1P1++];
        queueNumNow1--;
        
#ifdef USE_MPI_IO
        if (!cmd_info_->splitWrite_ && !out_is_zip_) {
            // MPI single file write mode
            long long my_offset;
            long long chunk_size = now.second;
            
            MPI_Fetch_and_op(&chunk_size, &my_offset, MPI_LONG_LONG, 
                            0, 0, MPI_SUM, mpi_win1_);
            MPI_Win_flush(0, mpi_win1_);
            
            MPI_Status status;
            MPI_File_write_at(mpi_file1_, my_offset, now.first, now.second, 
                             MPI_CHAR, &status);
            
            delete[] now.first;
        } else
#endif
        {
            if (out_is_zip_) {
                if (cmd_info_->use_pigz_) {
                    while (pigzQueueNumNow1 > pigzQueueSizeLim1) {
                        //debug_printf("waiting to push a chunk to pigz queue1\n");
                        usleep(100);
                    }
                    pigzQueue1->enqueue(now);
                    pigzQueueNumNow1++;
                } else {
                    int written = gzwrite(zip_out_stream1, now.first, now.second);
                    if (written != now.second) {
                        fprintf(stderr, "gzwrite error\n");
                        exit(0);
                    }
                    delete[] now.first;
                }
            } else {
                out_stream1_.write(now.first, now.second);
                delete[] now.first;
            }
        }
    }
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {

        } else {
            if (zip_out_stream1) {
                gzflush(zip_out_stream1, Z_FINISH);
                gzclose(zip_out_stream1);
                zip_out_stream1 = NULL;
            }
        }

    } else {
#ifdef USE_MPI_IO
        if (!cmd_info_->splitWrite_ && mpi_file1_ != MPI_FILE_NULL) {
            MPI_File_close(&mpi_file1_);
            MPI_Win_unlock_all(mpi_win1_);
            MPI_Win_free(&mpi_win1_);
            debug_printf("Rank %d: MPI file1 closed\n", my_rank_);
        } else
#endif
        {
            out_stream1_.close();
        }
    }
    debug_printf("write1 cost %.5f\n", GetTime() - t0);
}

/**
 * @brief a function to write pe data from out_data2 queue to file2
 */
void PeQc::WriteSeFastqTask2() {
    debug_printf("WriteSeFastqTask2 started\n");
    double t0 = GetTime();
    int cnt = 0;
    bool overWhile = 0;
    pair<char *, int> now;
    while (true) {
        while (queueNumNow2 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                break;
            }
            usleep(100);
        }
        if (overWhile) break;
        now = out_queue2_[queue2P1++];
        queueNumNow2--;
        
#ifdef USE_MPI_IO
        if (!cmd_info_->splitWrite_ && !out_is_zip_) {
            // MPI single file write mode
            long long my_offset;
            long long chunk_size = now.second;
            
            MPI_Fetch_and_op(&chunk_size, &my_offset, MPI_LONG_LONG, 
                            0, 0, MPI_SUM, mpi_win2_);
            MPI_Win_flush(0, mpi_win2_);
            
            MPI_Status status;
            MPI_File_write_at(mpi_file2_, my_offset, now.first, now.second, 
                             MPI_CHAR, &status);
            
            delete[] now.first;
        } else
#endif
        {
            if (out_is_zip_) {
                if (cmd_info_->use_pigz_) {
                    while (pigzQueueNumNow2 > pigzQueueSizeLim2) {
                        //debug_printf("waiting to push a chunk to pigz queue2\n");
                        usleep(100);
                    }
                    pigzQueue2->enqueue(now);
                    pigzQueueNumNow2++;
                } else {
                    int written = gzwrite(zip_out_stream2, now.first, now.second);
                    if (written != now.second) {
                        fprintf(stderr, "GG");
                        exit(0);
                    }

                    delete[] now.first;
                }
            } else {
                out_stream2_.write(now.first, now.second);
                delete[] now.first;
            }
        }
    }

    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {

        } else {
            if (zip_out_stream2) {
                gzflush(zip_out_stream2, Z_FINISH);
                gzclose(zip_out_stream2);
                zip_out_stream2 = NULL;
            }
        }

    } else {
#ifdef USE_MPI_IO
        if (!cmd_info_->splitWrite_ && mpi_file2_ != MPI_FILE_NULL) {
            MPI_File_close(&mpi_file2_);
            MPI_Win_unlock_all(mpi_win2_);
            MPI_Win_free(&mpi_win2_);
            debug_printf("Rank %d: MPI file2 closed\n", my_rank_);
        } else
#endif
        {
            out_stream2_.close();
        }
    }
    debug_printf("write2 cost %.5f\n", GetTime() - t0);
}

/*
void PeQc::PugzTask1() {
    debug_printf("pugz1 start\n");
    double t0 = GetTime();
    main_pugz(cmd_info_->in_file_name1_, cmd_info_->pugz_threads_, pugzQueue1, &producerDone);
    pugzDone1 = 1;
    debug_printf("pugz1 done, cost %.6f\n", GetTime() - t0);
}

void PeQc::PugzTask2() {
    debug_printf("pugz2 start\n");
    double t0 = GetTime();
    main_pugz(cmd_info_->in_file_name2_, cmd_info_->pugz_threads_, pugzQueue2, &producerDone);
    pugzDone2 = 1;
    debug_printf("pugz2 done, cost %.6f\n", GetTime() - t0);
}
*/

void PeQc::PugzTask1() {
    debug_printf("pragzip1 start\n");
    double t0 = GetTime();
 
    int cnt = 6;

    char **infos = new char *[6];
    infos[0] = "./pragzip";
    infos[1] = "-c";
    infos[2] = "-d";
    infos[3] = "-P";
    int th_num = cmd_info_->pugz_threads_;
    string th_num_s = to_string(th_num);
    infos[4] = new char[th_num_s.length() + 1];
    memcpy(infos[4], th_num_s.c_str(), th_num_s.length());
    infos[4][th_num_s.length()] = '\0';
    string in_file = cmd_info_->in_file_name1_;
    infos[5] = new char[in_file.length() + 1];
    memcpy(infos[5], in_file.c_str(), in_file.length());
    infos[5][in_file.length()] = '\0';

    main_pragzip(cnt, infos, pugzQueue1, &producerDone);

    pugzDone1 = 1;
    debug_printf("pragzip1 done, cost %.6f\n", GetTime() - t0);
}

void PeQc::PugzTask2() {
    debug_printf("pragzip2 start\n");
    double t0 = GetTime();

    int cnt = 6;

    char **infos = new char *[6];
    infos[0] = "./pragzip";
    infos[1] = "-c";
    infos[2] = "-d";
    infos[3] = "-P";
    int th_num = cmd_info_->pugz_threads_;
    string th_num_s = to_string(th_num);
    infos[4] = new char[th_num_s.length() + 1];
    memcpy(infos[4], th_num_s.c_str(), th_num_s.length());
    infos[4][th_num_s.length()] = '\0';
    string in_file = cmd_info_->in_file_name2_;
    infos[5] = new char[in_file.length() + 1];
    memcpy(infos[5], in_file.c_str(), in_file.length());
    infos[5][in_file.length()] = '\0';

    main_pragzip(cnt, infos, pugzQueue2, &producerDone);

    pugzDone2 = 1;
    debug_printf("pragzip2 done, cost %.6f\n", GetTime() - t0);
}


void PeQc::PigzTask1() {
    int cnt = 10;

    char **infos = new char *[10];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = cmd_info_->pigz_threads_;
    //    fprintf(stderr, "th num is %d\n", th_num);
    string th_num_s = to_string(th_num);
    //    fprintf(stderr, "th num s is %s\n", th_num_s.c_str());
    //    fprintf(stderr, "th num s len is %d\n", th_num_s.length());

    infos[2] = new char[th_num_s.length() + 1];
    memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";

    string tmp_level = to_string(cmd_info_->compression_level_);
    tmp_level = "-" + tmp_level;
    infos[4] = new char[tmp_level.length() + 1];
    memcpy(infos[4], tmp_level.c_str(), tmp_level.length());
    infos[4][tmp_level.length()] = '\0';


    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    string out_name1 = cmd_info_->out_file_name1_;
    string out_file = out_name1.substr(0, out_name1.find(".gz"));
    //    fprintf(stderr, "th out_file is %s\n", out_file.c_str());
    //    fprintf(stderr, "th out_file len is %d\n", out_file.length());
    infos[8] = new char[out_file.length() + 1];
    memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';
    infos[9] = "-v";
    main_pigz(cnt, infos, pigzQueue1, &writerDone1, pigzLast1, &pigzQueueNumNow1);
    debug_printf("pigz1 done\n");
}

void PeQc::PigzTask2() {

    int cnt = 10;

    char **infos = new char *[10];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = cmd_info_->pigz_threads_;
    //    fprintf(stderr, "th num is %d\n", th_num);
    string th_num_s = to_string(th_num);
    //    fprintf(stderr, "th num s is %s\n", th_num_s.c_str());
    //    fprintf(stderr, "th num s len is %d\n", th_num_s.length());

    infos[2] = new char[th_num_s.length() + 1];
    memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";

    string tmp_level = to_string(cmd_info_->compression_level_);
    tmp_level = "-" + tmp_level;
    infos[4] = new char[tmp_level.length() + 1];
    memcpy(infos[4], tmp_level.c_str(), tmp_level.length());
    infos[4][tmp_level.length()] = '\0';

    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    string out_name2 = cmd_info_->out_file_name2_;
    string out_file = out_name2.substr(0, out_name2.find(".gz"));
    //    fprintf(stderr, "th out_file is %s\n", out_file.c_str());
    //    fprintf(stderr, "th out_file len is %d\n", out_file.length());
    infos[8] = new char[out_file.length() + 1];
    memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';
    infos[9] = "-v";
    main_pigz(cnt, infos, pigzQueue2, &writerDone2, pigzLast2, &pigzQueueNumNow2);
    debug_printf("pigz2 done\n");
}

/**
 * @brief do QC for pair-end data
 */
void PeQc::ProcessPeFastq() {
    if (cmd_info_->interleaved_in_) {
        auto *fastqPool = new rabbit::fq::FastqDataPool(64, 1 << 22);
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);
        auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            p_thread_info[t] = new ThreadInfo(cmd_info_, true);
        }
        thread *write_thread1;
        thread *write_thread2;
        if (cmd_info_->write_data_) {
            write_thread1 = new thread(bind(&PeQc::WriteSeFastqTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2 = new thread(bind(&PeQc::WriteSeFastqTask2, this));
        }

        thread producer(
                bind(&PeQc::ProducerPeInterFastqTask, this, cmd_info_->in_file_name1_, fastqPool,
                    ref(queue1)));
        auto **threads = new thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new thread(
                    bind(&PeQc::ConsumerPeInterFastqTask, this, p_thread_info[t], fastqPool, ref(queue1)));
        }
        producer.join();
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t]->join();
        }
        if (cmd_info_->write_data_) {
            write_thread1->join();
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2->join();
        }
        debug_printf("all thread done\n");
        debug_printf("now merge thread info\n");
        vector<State *> pre_vec_state1;
        vector<State *> pre_vec_state2;
        vector<State *> aft_vec_state1;
        vector<State *> aft_vec_state2;

        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            pre_vec_state1.push_back(p_thread_info[t]->pre_state1_);
            pre_vec_state2.push_back(p_thread_info[t]->pre_state2_);
            aft_vec_state1.push_back(p_thread_info[t]->aft_state1_);
            aft_vec_state2.push_back(p_thread_info[t]->aft_state2_);
        }
        auto pre_state1_tmp = State::MergeStates(pre_vec_state1);
        auto pre_state2_tmp = State::MergeStates(pre_vec_state2);
        auto aft_state1_tmp = State::MergeStates(aft_vec_state1);
        auto aft_state2_tmp = State::MergeStates(aft_vec_state2);

        debug_printf("merge done\n");

#ifdef USE_MPI_IO
        // MPI: merge state information across all processes
        vector<State *> pre_state1_mpis, pre_state2_mpis, aft_state1_mpis, aft_state2_mpis;
        pre_state1_mpis.push_back(pre_state1_tmp);
        pre_state2_mpis.push_back(pre_state2_tmp);
        aft_state1_mpis.push_back(aft_state1_tmp);
        aft_state2_mpis.push_back(aft_state2_tmp);

        // Gather pre_state1
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
                delete[] info;
                pre_state1_mpis.push_back(tmp_state);
            }
        } else {
            string pre_state1_is = pre_state1_tmp->ParseString();
            int now_size = pre_state1_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(pre_state1_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gather pre_state2
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, true);
                delete[] info;
                pre_state2_mpis.push_back(tmp_state);
            }
        } else {
            string pre_state2_is = pre_state2_tmp->ParseString();
            int now_size = pre_state2_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(pre_state2_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gather aft_state1
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
                delete[] info;
                aft_state1_mpis.push_back(tmp_state);
            }
        } else {
            string aft_state1_is = aft_state1_tmp->ParseString();
            int now_size = aft_state1_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(aft_state1_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gather aft_state2
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, true);
                delete[] info;
                aft_state2_mpis.push_back(tmp_state);
            }
        } else {
            string aft_state2_is = aft_state2_tmp->ParseString();
            int now_size = aft_state2_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(aft_state2_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Non-zero ranks clean up and exit
        if(my_rank_ != 0) {
            delete pre_state1_tmp;
            delete pre_state2_tmp;
            delete aft_state1_tmp;
            delete aft_state2_tmp;
            delete fastqPool;
            for (int t = 0; t < cmd_info_->thread_number_; t++) {
                delete p_thread_info[t];
            }
            delete[] p_thread_info;
            return;
        }

        // Only rank 0 continues: merge all states
        auto pre_state1 = State::MergeStates(pre_state1_mpis);
        auto pre_state2 = State::MergeStates(pre_state2_mpis);
        auto aft_state1 = State::MergeStates(aft_state1_mpis);
        auto aft_state2 = State::MergeStates(aft_state2_mpis);
        for(int i = 1; i < comm_size_; i++) {
            delete pre_state1_mpis[i];
            delete pre_state2_mpis[i];
            delete aft_state1_mpis[i];
            delete aft_state2_mpis[i];
        }
        delete pre_state1_tmp;
        delete pre_state2_tmp;
        delete aft_state1_tmp;
        delete aft_state2_tmp;
#else
        auto pre_state1 = pre_state1_tmp;
        auto pre_state2 = pre_state2_tmp;
        auto aft_state1 = aft_state1_tmp;
        auto aft_state2 = aft_state2_tmp;
#endif

        fprintf(stderr, "\nprint read1 (before filter) info :\n");
        State::PrintStates(pre_state1);
        fprintf(stderr, "\nprint read1 (after filter) info :\n");
        State::PrintStates(aft_state1);
        fprintf(stderr, "\n");

        fprintf(stderr, "\nprint read2 (before filter) info :\n");
        State::PrintStates(pre_state2);
        fprintf(stderr, "\nprint read2 (after filter) info :\n");
        State::PrintStates(aft_state2);
        fprintf(stderr, "\n");
        if (cmd_info_->print_what_trimmed_) {
            State::PrintAdapterToFile(aft_state1);
            State::PrintAdapterToFile(aft_state2);
        }
        State::PrintFilterResults(aft_state1);
        fprintf(stderr, "\n");

        if (cmd_info_->do_overrepresentation_ && cmd_info_->print_ORP_seqs_) {

            auto pre_hash_graph1 = pre_state1->GetHashGraph();
            int pre_hash_num1 = pre_state1->GetHashNum();
            auto pre_hash_graph2 = pre_state2->GetHashGraph();
            int pre_hash_num2 = pre_state2->GetHashNum();

            auto aft_hash_graph1 = aft_state1->GetHashGraph();
            int aft_hash_num1 = aft_state1->GetHashNum();
            auto aft_hash_graph2 = aft_state2->GetHashGraph();
            int aft_hash_num2 = aft_state2->GetHashNum();


            int spg = cmd_info_->overrepresentation_sampling_;
            ofstream ofs;

            string srr_name1 = cmd_info_->in_file_name1_;
            srr_name1 = PaseFileName(srr_name1);

            string srr_name2 = cmd_info_->in_file_name2_;
            srr_name2 = PaseFileName(srr_name2);

            string out_name1 = "pe_" + srr_name1 + "_before_ORP_sequences.txt";
            ofs.open(out_name1, ifstream::out);

            int cnt1 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < pre_hash_num1; i++) {
                if (!overRepPassed(pre_hash_graph1[i].seq, pre_hash_graph1[i].cnt, spg)) continue;
                ofs << pre_hash_graph1[i].seq << " " << pre_hash_graph1[i].cnt << "\n";
                cnt1++;
            }
            ofs.close();
            fprintf(stderr, "in %s (before filter) find %d possible overrepresented sequences (store in %s)\n",
                    srr_name1.c_str(), cnt1, out_name1.c_str());

            string out_name2 = "pe_" + srr_name2 + "_before_ORP_sequences.txt";
            ofs.open(out_name2, ifstream::out);
            int cnt2 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < pre_hash_num2; i++) {
                if (!overRepPassed(pre_hash_graph2[i].seq, pre_hash_graph2[i].cnt, spg)) continue;
                ofs << pre_hash_graph2[i].seq << " " << pre_hash_graph2[i].cnt << "\n";
                cnt2++;
            }
            ofs.close();
            fprintf(stderr, "in %s (before filter) find %d possible overrepresented sequences (store in %s)\n",
                    srr_name2.c_str(), cnt2, out_name2.c_str());


            out_name1 = "pe_" + srr_name1 + "_after_ORP_sequences.txt";
            ofs.open(out_name1, ifstream::out);
            cnt1 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < aft_hash_num1; i++) {
                if (!overRepPassed(aft_hash_graph1[i].seq, aft_hash_graph1[i].cnt, spg)) continue;
                ofs << aft_hash_graph1[i].seq << " " << aft_hash_graph1[i].cnt << "\n";
                cnt1++;
            }
            ofs.close();
            fprintf(stderr, "in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name1.c_str(),
                    cnt1, out_name1.c_str());

            out_name2 = "pe_" + srr_name2 + "_after_ORP_sequences.txt";
            ofs.open(out_name2, ifstream::out);
            cnt2 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < aft_hash_num2; i++) {
                if (!overRepPassed(aft_hash_graph2[i].seq, aft_hash_graph2[i].cnt, spg)) continue;
                ofs << aft_hash_graph2[i].seq << " " << aft_hash_graph2[i].cnt << "\n";
                cnt2++;
            }
            ofs.close();
            fprintf(stderr, "in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name2.c_str(),
                    cnt2, out_name2.c_str());

            fprintf(stderr, "\n");
        }
        int *dupHist = NULL;
        double *dupMeanGC = NULL;
        double dupRate = 0.0;
        int histSize = 32;
        if (cmd_info_->state_duplicate_) {
            dupHist = new int[histSize];
            memset(dupHist, 0, sizeof(int) * histSize);
            dupMeanGC = new double[histSize];
            memset(dupMeanGC, 0, sizeof(double) * histSize);
            dupRate = duplicate_->statAll(dupHist, dupMeanGC, histSize);
            fprintf(stderr, "Duplication rate : %.5f %%\n", dupRate * 100.0);
            delete[] dupHist;
            delete[] dupMeanGC;
        }

        int64_t *merge_insert_size;
        if (!cmd_info_->no_insert_size_) {
            merge_insert_size = new int64_t[cmd_info_->max_insert_size_ + 1];
            memset(merge_insert_size, 0, sizeof(int64_t) * (cmd_info_->max_insert_size_ + 1));

            for (int t = 0; t < cmd_info_->thread_number_; t++) {
                for (int i = 0; i <= cmd_info_->max_insert_size_; i++) {
                    merge_insert_size[i] += p_thread_info[t]->insert_size_dist_[i];
                }
            }
            int mx_id = 0;
            for (int i = 0; i < cmd_info_->max_insert_size_; i++) {
                if (merge_insert_size[i] > merge_insert_size[mx_id]) mx_id = i;
            }
            //fprintf(stderr, "Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
            fprintf(stderr, "Insert size peak (based on PE overlap analyze): %d\n", mx_id);
        }
        string srr_name1 = cmd_info_->in_file_name1_;
        srr_name1 = PaseFileName(srr_name1);
        string srr_name2 = cmd_info_->in_file_name2_;
        srr_name2 = PaseFileName(srr_name2);
        Repoter::ReportHtmlPe(srr_name1 + "_" + srr_name2 + "_RabbitQCPlus.html", pre_state1, pre_state2, aft_state1,
                aft_state2, cmd_info_->in_file_name1_,
                cmd_info_->in_file_name2_, dupRate * 100.0, merge_insert_size);
        debug_printf("report done\n");

        delete pre_state1;
        delete pre_state2;
        delete aft_state1;
        delete aft_state2;
        if (!cmd_info_->no_insert_size_)
            delete[] merge_insert_size;


        delete fastqPool;
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            delete p_thread_info[t];
            delete threads[t];
        }
        delete[] threads;
        delete[] p_thread_info;
        if (cmd_info_->write_data_) {
            delete write_thread1;
            if (cmd_info_->interleaved_out_ == 0)
                delete write_thread2;
        }
    } else {
        thread *carer;
        if(cmd_info_->do_correction_with_care_) {
            carer = new thread(bind(&PeQc::careProcess, this));
            //cmd_info_->in_file_name1_ = "./tmp.fq";
            while(careStartWrite == 0) {
                usleep(10000);
            }
            fprintf(stderr, "QC start...\n");
            if(changeNum == 0) {
                carer->join();
                delete carer;
                cmd_info_->do_correction_with_care_ = 0;
            } else {
                cmd_info_->use_pugz_ = 0;
            }
        }
        if(cmd_info_->use_pigz_) {
            if(cmd_info_->do_correction_with_care_) {
                carer->join();
                delete carer;
            }
        }


        thread *pugzer1;
        thread *pugzer2;

        if (cmd_info_->use_pugz_) {
            pugzer1 = new thread(bind(&::PeQc::PugzTask1, this));
            pugzer2 = new thread(bind(&::PeQc::PugzTask2, this));
        }


        auto *fastqPool = new rabbit::fq::FastqDataPool(64, 1 << 22);
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(64, 1);
        auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            p_thread_info[t] = new ThreadInfo(cmd_info_, true);
        }
        thread *write_thread1;
        thread *write_thread2;
        if (cmd_info_->write_data_) {
            write_thread1 = new thread(bind(&PeQc::WriteSeFastqTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2 = new thread(bind(&PeQc::WriteSeFastqTask2, this));
        }
        thread *pigzer1;
        thread *pigzer2;
        if (cmd_info_->use_pigz_) {
            pigzer1 = new thread(bind(&PeQc::PigzTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                pigzer2 = new thread(bind(&PeQc::PigzTask2, this));
        }
        thread producer(
                bind(&PeQc::ProducerPeFastqTask, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_,
                    fastqPool, ref(queue1)));
        auto **threads = new thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new thread(
                    bind(&PeQc::ConsumerPeFastqTask, this, p_thread_info[t], fastqPool, ref(queue1)));
        }

        if (cmd_info_->use_pugz_) {
            pugzer1->join();
            pugzer2->join();
        }
        if(cmd_info_->use_pigz_ == 0) {
            if(cmd_info_->do_correction_with_care_) {
                carer->join();
                delete carer;
            }
        }
        

        producer.join();
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t]->join();
        }

        if (cmd_info_->write_data_) {
            write_thread1->join();
            writerDone1 = 1;
            if (cmd_info_->interleaved_out_ == 0) {

                write_thread2->join();
                writerDone2 = 1;
            }
        }

        if (cmd_info_->use_pigz_) {
            pigzer1->join();
            if (cmd_info_->interleaved_out_ == 0)
                pigzer2->join();
        }
        debug_printf("all thread done\n");
        debug_printf("now merge thread info\n");
        vector<State *> pre_vec_state1;
        vector<State *> pre_vec_state2;
        vector<State *> aft_vec_state1;
        vector<State *> aft_vec_state2;

        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            pre_vec_state1.push_back(p_thread_info[t]->pre_state1_);
            pre_vec_state2.push_back(p_thread_info[t]->pre_state2_);
            aft_vec_state1.push_back(p_thread_info[t]->aft_state1_);
            aft_vec_state2.push_back(p_thread_info[t]->aft_state2_);
        }
        auto pre_state1_tmp = State::MergeStates(pre_vec_state1);
        auto pre_state2_tmp = State::MergeStates(pre_vec_state2);
        auto aft_state1_tmp = State::MergeStates(aft_vec_state1);
        auto aft_state2_tmp = State::MergeStates(aft_vec_state2);

#ifdef USE_MPI_IO
        // MPI: merge state information across all processes
        vector<State *> pre_state1_mpis, pre_state2_mpis, aft_state1_mpis, aft_state2_mpis;
        pre_state1_mpis.push_back(pre_state1_tmp);
        pre_state2_mpis.push_back(pre_state2_tmp);
        aft_state1_mpis.push_back(aft_state1_tmp);
        aft_state2_mpis.push_back(aft_state2_tmp);

        // Gather pre_state1
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
                delete[] info;
                pre_state1_mpis.push_back(tmp_state);
            }
        } else {
            string pre_state1_is = pre_state1_tmp->ParseString();
            int now_size = pre_state1_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(pre_state1_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gather pre_state2
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, true);
                delete[] info;
                pre_state2_mpis.push_back(tmp_state);
            }
        } else {
            string pre_state2_is = pre_state2_tmp->ParseString();
            int now_size = pre_state2_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(pre_state2_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gather aft_state1
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
                delete[] info;
                aft_state1_mpis.push_back(tmp_state);
            }
        } else {
            string aft_state1_is = aft_state1_tmp->ParseString();
            int now_size = aft_state1_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(aft_state1_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gather aft_state2
        if(my_rank_ == 0) {
            for(int i = 1; i < comm_size_; i++) {
                int now_size = 0;
                MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                char* info = new char[now_size];
                MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, true);
                delete[] info;
                aft_state2_mpis.push_back(tmp_state);
            }
        } else {
            string aft_state2_is = aft_state2_tmp->ParseString();
            int now_size = aft_state2_is.length();
            MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(aft_state2_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Non-zero ranks clean up and exit
        if(my_rank_ != 0) {
            delete pre_state1_tmp;
            delete pre_state2_tmp;
            delete aft_state1_tmp;
            delete aft_state2_tmp;
            delete fastqPool;
            for (int t = 0; t < cmd_info_->thread_number_; t++) {
                delete threads[t];
                delete p_thread_info[t];
            }
            delete[] threads;
            delete[] p_thread_info;
            if (cmd_info_->write_data_) {
                delete write_thread1;
                if (cmd_info_->interleaved_out_ == 0)
                    delete write_thread2;
            }
            return;
        }

        // Only rank 0 continues: merge all states
        auto pre_state1 = State::MergeStates(pre_state1_mpis);
        auto pre_state2 = State::MergeStates(pre_state2_mpis);
        auto aft_state1 = State::MergeStates(aft_state1_mpis);
        auto aft_state2 = State::MergeStates(aft_state2_mpis);
        for(int i = 1; i < comm_size_; i++) {
            delete pre_state1_mpis[i];
            delete pre_state2_mpis[i];
            delete aft_state1_mpis[i];
            delete aft_state2_mpis[i];
        }
        delete pre_state1_tmp;
        delete pre_state2_tmp;
        delete aft_state1_tmp;
        delete aft_state2_tmp;
#else
        auto pre_state1 = pre_state1_tmp;
        auto pre_state2 = pre_state2_tmp;
        auto aft_state1 = aft_state1_tmp;
        auto aft_state2 = aft_state2_tmp;
#endif

        if (cmd_info_->do_overrepresentation_) {
            debug_printf("orp cost %f\n", pre_state1->GetOrpCost() + pre_state2->GetOrpCost() + aft_state1->GetOrpCost() + aft_state2->GetOrpCost());
        }
        fprintf(stderr, "merge done\n");
        fprintf(stderr, "\nprint read1 (before filter) info :\n");
        State::PrintStates(pre_state1);
        fprintf(stderr, "\nprint read1 (after filter) info :\n");
        State::PrintStates(aft_state1);
        fprintf(stderr, "\n");

        fprintf(stderr, "\nprint read2 (before filter) info :\n");
        State::PrintStates(pre_state2);
        fprintf(stderr, "\nprint read2 (after filter) info :\n");
        State::PrintStates(aft_state2);
        fprintf(stderr, "\n");
        if (cmd_info_->print_what_trimmed_) {
            State::PrintAdapterToFile(aft_state1);
            State::PrintAdapterToFile(aft_state2);
        }
        State::PrintFilterResults(aft_state1);
        fprintf(stderr, "\n");
        if (cmd_info_->do_overrepresentation_ && cmd_info_->print_ORP_seqs_) {

            auto pre_hash_graph1 = pre_state1->GetHashGraph();
            int pre_hash_num1 = pre_state1->GetHashNum();
            auto pre_hash_graph2 = pre_state2->GetHashGraph();
            int pre_hash_num2 = pre_state2->GetHashNum();

            auto aft_hash_graph1 = aft_state1->GetHashGraph();
            int aft_hash_num1 = aft_state1->GetHashNum();
            auto aft_hash_graph2 = aft_state2->GetHashGraph();
            int aft_hash_num2 = aft_state2->GetHashNum();


            int spg = cmd_info_->overrepresentation_sampling_;
            ofstream ofs;

            string srr_name1 = cmd_info_->in_file_name1_;
            srr_name1 = PaseFileName(srr_name1);

            string srr_name2 = cmd_info_->in_file_name2_;
            srr_name2 = PaseFileName(srr_name2);

            string out_name1 = "pe_" + srr_name1 + "_before_ORP_sequences.txt";
            ofs.open(out_name1, ifstream::out);

            int cnt1 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < pre_hash_num1; i++) {
                if (!overRepPassed(pre_hash_graph1[i].seq, pre_hash_graph1[i].cnt, spg)) continue;
                ofs << pre_hash_graph1[i].seq << " " << pre_hash_graph1[i].cnt << "\n";
                cnt1++;
            }
            ofs.close();
            fprintf(stderr, "in %s (before filter) find %d possible overrepresented sequences (store in %s)\n",
                    srr_name1.c_str(), cnt1, out_name1.c_str());

            string out_name2 = "pe_" + srr_name2 + "_before_ORP_sequences.txt";
            ofs.open(out_name2, ifstream::out);
            int cnt2 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < pre_hash_num2; i++) {
                if (!overRepPassed(pre_hash_graph2[i].seq, pre_hash_graph2[i].cnt, spg)) continue;
                ofs << pre_hash_graph2[i].seq << " " << pre_hash_graph2[i].cnt << "\n";
                cnt2++;
            }
            ofs.close();
            fprintf(stderr, "in %s (before filter) find %d possible overrepresented sequences (store in %s)\n",
                    srr_name2.c_str(), cnt2, out_name2.c_str());


            out_name1 = "pe_" + srr_name1 + "_after_ORP_sequences.txt";
            ofs.open(out_name1, ifstream::out);
            cnt1 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < aft_hash_num1; i++) {
                if (!overRepPassed(aft_hash_graph1[i].seq, aft_hash_graph1[i].cnt, spg)) continue;
                ofs << aft_hash_graph1[i].seq << " " << aft_hash_graph1[i].cnt << "\n";
                cnt1++;
            }
            ofs.close();
            fprintf(stderr, "in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name1.c_str(),
                    cnt1, out_name1.c_str());

            out_name2 = "pe_" + srr_name2 + "_after_ORP_sequences.txt";
            ofs.open(out_name2, ifstream::out);
            cnt2 = 0;
            ofs << "sequence count"
                << "\n";
            for (int i = 0; i < aft_hash_num2; i++) {
                if (!overRepPassed(aft_hash_graph2[i].seq, aft_hash_graph2[i].cnt, spg)) continue;
                ofs << aft_hash_graph2[i].seq << " " << aft_hash_graph2[i].cnt << "\n";
                cnt2++;
            }
            ofs.close();
            fprintf(stderr, "in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name2.c_str(),
                    cnt2, out_name2.c_str());

            fprintf(stderr, "\n");
        }
        int *dupHist = NULL;
        double *dupMeanGC = NULL;
        double dupRate = 0.0;
        int histSize = 32;
        if (cmd_info_->state_duplicate_) {
            dupHist = new int[histSize];
            memset(dupHist, 0, sizeof(int) * histSize);
            dupMeanGC = new double[histSize];
            memset(dupMeanGC, 0, sizeof(double) * histSize);
            dupRate = duplicate_->statAll(dupHist, dupMeanGC, histSize);
            fprintf(stderr, "Duplication rate : %.5f %%\n", dupRate * 100.0);
            delete[] dupHist;
            delete[] dupMeanGC;
        }

        int64_t *merge_insert_size;
        if (!cmd_info_->no_insert_size_) {
            merge_insert_size = new int64_t[cmd_info_->max_insert_size_ + 1];
            memset(merge_insert_size, 0, sizeof(int64_t) * (cmd_info_->max_insert_size_ + 1));

            for (int t = 0; t < cmd_info_->thread_number_; t++) {
                for (int i = 0; i <= cmd_info_->max_insert_size_; i++) {
                    merge_insert_size[i] += p_thread_info[t]->insert_size_dist_[i];
                }
            }
            int mx_id = 0;
            for (int i = 0; i < cmd_info_->max_insert_size_; i++) {
                if (merge_insert_size[i] > merge_insert_size[mx_id]) mx_id = i;
            }
            //fprintf(stderr, "Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
            fprintf(stderr, "Insert size peak (based on PE overlap analyze): %d\n", mx_id);
        }
        string srr_name1 = cmd_info_->in_file_name1_;
        srr_name1 = PaseFileName(srr_name1);
        string srr_name2 = cmd_info_->in_file_name2_;
        srr_name2 = PaseFileName(srr_name2);
        Repoter::ReportHtmlPe(srr_name1 + "_" + srr_name2 + "_RabbitQCPlus.html", pre_state1, pre_state2, aft_state1,
                aft_state2, cmd_info_->in_file_name1_,
                cmd_info_->in_file_name2_, dupRate * 100.0, merge_insert_size);

        JsonReporter::ReportPe(srr_name1 + "_" + srr_name2 + "_RabbitQCPlus.json", pre_state1, pre_state2, aft_state1,
                               aft_state2, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_, dupRate * 100.0, merge_insert_size);
        debug_printf("report done\n");
        delete pre_state1;
        delete pre_state2;
        delete aft_state1;
        delete aft_state2;
        if (!cmd_info_->no_insert_size_)
            delete[] merge_insert_size;

        delete fastqPool;
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            delete p_thread_info[t];
            delete threads[t];
        }

        delete[] threads;
        delete[] p_thread_info;
        if (cmd_info_->write_data_) {
            delete write_thread1;
            if (cmd_info_->interleaved_out_ == 0)
                delete write_thread2;
        }

    }
}
