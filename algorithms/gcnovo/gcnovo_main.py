# Copyright 2019 Rui Qiao. All Rights Reserved.
#
# DeepNovoV2 is publicly available for non-commercial uses.
# ==============================================================================
import time
import torch
import logging
import logging.config
import config
import os

from mgf2feature import mgftofeature
from train_func import train, build_model, validation, perplexity
# from train_func2 import train, build_model, validation, perplexity
# from train_func import train, build_model, validation, perplexity
from data_reader import DeepNovoDenovoDataset, collate_func, DeepNovoTrainDataset
from model_gcn import InferenceModelWrapper
from denovo import IonCNNDenovo
# from model import InferenceModelWrapper
# from denovo import IonCNNDenovo
# from model_gcn2 import InferenceModelWrapper
# from denovo2 import IonCNNDenovo
from writer import DenovoWriter
from init_args import init_args
import worker_test
# from dia_script_select import find_score_cutoff
import datetime
from worker_test import WorkerTest
logger = logging.getLogger(__name__)

def engine_1(args):
    # train + search denovo + test
    print(f"training mode")
    torch.cuda.empty_cache()
    train(args=args)
    """
    search denovo
    """
    torch.cuda.empty_cache()
    start = time.time()
    print("denovo mode")
    data_reader = DeepNovoDenovoDataset(feature_filename=args.denovo_input_feature_file,
                                        spectrum_filename=args.denovo_input_spectrum_file,
                                        args=args)
    denovo_worker = IonCNNDenovo(args=args)
    # forward_deepnovo, backward_deepnovo, init_net = build_model(training=False)
    # model_wrapper = InferenceModelWrapper(forward_deepnovo, backward_deepnovo, init_net)
    forward_deepnovo, backward_deepnovo, init_net = build_model(args=args, training=False)
    model_wrapper = InferenceModelWrapper(forward_deepnovo, backward_deepnovo, init_net)
    writer = DenovoWriter(args=args)
    denovo_worker.search_denovo(model_wrapper, data_reader, writer)
    torch.cuda.empty_cache()
    print('using time:', time.time() - start)

    print("test mode")
    worker_test = WorkerTest(args=args)
    worker_test.test_accuracy()

    # show 95 accuracy score threshold
    accuracy_cutoff = 0.95
    accuracy_file = args.accuracy_file
    # score_cutoff = find_score_cutoff(accuracy_file, accuracy_cutoff)

def engine_2(args):
    # search denovo + test
    """
    search denovo
    """
    mgftofeature(args.denovo_input_spectrum_file)
    torch.cuda.empty_cache()
    start = time.time()
    print("denovo mode")
    data_reader = DeepNovoDenovoDataset(feature_filename=args.denovo_input_feature_file,
                                        spectrum_filename=args.denovo_input_spectrum_file,
                                        args=args)
    denovo_worker = IonCNNDenovo(args=args)
    # forward_deepnovo, backward_deepnovo, init_net = build_model(training=False)
    # model_wrapper = InferenceModelWrapper(forward_deepnovo, backward_deepnovo, init_net)
    forward_deepnovo, backward_deepnovo, init_net = build_model(args=args, training=False)
    model_wrapper = InferenceModelWrapper(forward_deepnovo, backward_deepnovo, init_net)
    writer = DenovoWriter(args=args)
    denovo_worker.search_denovo(model_wrapper, data_reader, writer)
    torch.cuda.empty_cache()
    print('using time:', time.time() - start)

    # print("test mode")
    # worker_test = worker_test.WorkerTest(args=args)
    # worker_test.test_accuracy()
    # 
    # # show 95 accuracy score threshold
    # accuracy_cutoff = 0.95
    # accuracy_file = args.accuracy_file
    # # score_cutoff = find_score_cutoff(accuracy_file, accuracy_cutoff)

def engine_3(args):
    # test
    print("test mode")
    worker_test = WorkerTest(args=args)
    worker_test.test_accuracy()

    # show 95 accuracy score threshold
    accuracy_cutoff = 0.95
    accuracy_file = args.accuracy_file
    # score_cutoff = find_score_cutoff(accuracy_file, accuracy_cutoff)

def init_log(log_file_name):
    d = {
        'version': 1,
        'disable_existing_loggers': False,  # this fixes the problem
        'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
            },
        },
        'handlers': {
            'console': {
                'level': 'INFO',
                'class': 'logging.StreamHandler',
                'formatter': 'standard',
            },
            'file': {
                'level': 'DEBUG',
                'class': 'logging.FileHandler',
                'filename': log_file_name,
                'mode': 'w',
                'formatter': 'standard',
            }
        },
        'root': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
        }
    }
    logging.config.dictConfig(d)

if __name__ == '__main__':
    param_path = "/algo/param/params.cfg"
    # log_path = "./GCNovo"
    print(param_path)
    if os.path.isfile(param_path):
        # log_path += (param_path.split("/")[-1] + "_")
        dir, param_file = os.path.split(param_path)
        # log_file_name = "top5_" + param_file[-4] + ".log"
        now = datetime.datetime.now().strftime("%Y%m%d%H%M")
        args = init_args(param_path)
        # log_file_name = "./log/" + now + "(" + str(args.engine_model) + ").log"
        # log_file_name = log_path + now + "(" + str(args.engine_model) + ").log"
        # init_log(log_file_name=log_file_name)
        if os.path.exists(args.train_dir):
            pass
        else:
            os.makedirs(args.train_dir)
        if args.engine_model == 1:
            # print("engine model 1")
            engine_1(args=args)
        elif args.engine_model == 2:
            engine_2(args=args)
            #执行修改output格式代码
            import output_mapper
            output_mapper.output(args.denovo_output_file)
        elif args.engine_model == 3:
            engine_3(args=args)
    elif os.path.isdir(param_path):
        list_dir = os.listdir(param_path)
        list_dir.sort(key=lambda x: int(x[33]))
        print(list_dir)
        for file in list_dir:
            one_param_path = os.path.join(param_path, file)
            if os.path.isfile(one_param_path):
                now = datetime.datetime.now().strftime("%Y%m%d%H%M")
                args = init_args(one_param_path)
                # log_file_name = log_path + file+"_" + now + "(" + str(args.engine_model) + ").log"
                # init_log(log_file_name=log_file_name)
                if os.path.exists(args.train_dir):
                    pass
                else:
                    os.makedirs(args.train_dir)
                if args.engine_model == 1:
                    engine_1(args=args)
                elif args.engine_model == 2:
                    engine_2(args=args)
                elif args.engine_model == 3:
                    engine_3(args=args)
