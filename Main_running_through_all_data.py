import Main as M
import Cell_line_eater
import os
import time
import read_in_dan_descriptor_make_sam_from_graphs as Dan_fromat

# Dan functions (copy of crossvalidate classification)
import sys
from IPython import embed
import numpy as np
from scipy import interp
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc, f1_score
from sklearn.cross_validation import StratifiedKFold, StratifiedShuffleSplit
from sklearn import tree, ensemble, naive_bayes, neural_network, svm
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score, roc_curve
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.AtomPairs import Pairs, Torsions
from multiprocessing import Pool
import multiprocessing
import math as m
import csv
from IPython import embed
from operator import itemgetter
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn.metrics import accuracy_score, precision_recall_curve, matthews_corrcoef
from sklearn.metrics import average_precision_score, precision_score, recall_score
from sklearn import preprocessing
from sklearn import cross_validation
import scipy as scp
from functions import getXy, continuous_tanimoto, findScoreVersusAbsoluteAccuracy, featureRank
from scipy.stats import ranksums, spearmanr
import cPickle
from sklearn import metrics
from sklearn import ensemble, svm, preprocessing, feature_selection
from matplotlib import pyplot as p
import math as m
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.grid_search import RandomizedSearchCV, GridSearchCV

#Globals:
graph = "ominpath.ncol"#used graph
folder = "/home/dm729/ucc-fileserver/PycharmProjects/Graph_similarity_measures/" #used folder
sample_file = "dan_combination_train_signalink_targets_and_neighbours.csv" # To make accoruate file format using this file
Krishna_file = "cell_line_tissue_dictionary.txt" #Tissue translating if something is missing than this is
#corrigated based here the tissue annotation. Curently after the graph based influence done.

#Import annotations
SwissProtUniProtIds = M.uniprotin(folder+"uniprot-homo+sapiens.tab")
gene_name_uniprot_library = M.chip_annotation_1_to_chip_annotation_2(folder+"GPL13667-15572_annotation.csv", 14, 21, SwissProtUniProtIds)
print "Annotations are ready they are in the memory"

#First line creats the expressions of the cell lines. Theese are the so called .expr files.
# Theese are the list of expressing mRNAs in the cell lines. The graph vbased analyisis is always this is the first
# step. The function have two type of argument  1:
# SD as -1 SD expression do not counted as expressed
# PT as Percentage of expression less than 25 counted as not expressed.
# This will be first SD based and than percentage based.
Cell_line_eater.cell_line_extractor(folder, "gex.csv", "SD")
print "Cell lines data are created. "
expressions=[]
for each in os.listdir(folder):
    if each.endswith("SD.expr"):
        expressions += [each]

def run_through_expression(neigbourhood, propagation_type, expressions, expression_type):
    """
    This function creates the cell line graph based descriptions according to the paramters
    :param neigbourhood: int, tells how far should igraph chek the neighbours
    :param propagation_type: 1 if a vertex is reached from an oter one added to used vertexes
    2 add to used vertexes only the end of the cycle to prevent infomration backflow.
    3 not add to used vertexes, allow backword infomration flow
    :param expressions: expression files containing list
    :return:
    """
    for each in expressions:
        M.graph_from_expression_file_graph(folder+each, graph, gene_name_uniprot_library, neigbourhood, propagation_type,
                                           expression_type)
        b = time.clock()
        print "Cell line completed:", each, "Time ellapsed since start:", (b-a)/60, "minutes"
        print each, "done"
    gene_expression_descriptor_end_tag = "cell_line_gene_distance_affy_translation_only_SP_" + expression_type + ".celist"
    final_file_name = graph.replace(".ncol","_") +"type"+str(propagation_type)+"_neighborhood"+str(neigbourhood)+"_expressiontype"
    +expression_type+".csv"
    Dan_fromat.construct_descriptors_for_additional_cell_lines(folder+Krishna_file, "\t", folder, gene_expression_descriptor_end_tag)
    dan_descriptor_file = Dan_fromat.read_in_dan_descriptor(folder+sample_file)
    Dan_fromat.open_cell_lines_write_out_descriptor(folder, dan_descriptor_file, gene_expression_descriptor_end_tag,
                                     final_file_name)
    #This dictionarry contains the fil names which are in the modell.
    out = open("output_log.txt", "wb")
    fileNames = {
    #"Graph_descritpros" : input_dir+"dezso_combination_test_Signor_2nd_neighbors_GeneExpression_Mean_minus_one_SD_not_exp_mean.csv",
    #"0_additional_descriptors" : input_dir+prefix+'combination_'+dataset_type+add_in+'_additional_descriptors.csv',
    #"annotated_cell_line" : input_dir+prefix+'combination_'+dataset_type+add_in+'_annotated_cell_line.csv',
    #"annotated_disease" : input_dir+prefix+'combination_'+dataset_type+add_in+'_annotated_disease.csv',
    "annotated_targets" : input_dir+prefix+'combination_'+dataset_type+add_in+'_annotated_targets.csv',
    #"annotated_tissue" : input_dir+prefix+'combination_'+dataset_type+add_in+'_annotated_tissue.csv',
    #"predicted_targets" : input_dir+prefix+'combination_'+dataset_type+add_in+'_predicted_targets.csv',
    #"signalink_pathways" : input_dir+prefix+'combination_'+dataset_type+add_in+'_signalink_pathways.csv',
    #"signalink_targets_and_neighbours" : prefix+'combination_'+dataset_type+add_in+'_signalink_targets_and_neighbours.csv',
    #"signalink_targets_and_neighbours" : 'neighbours_and_targets.csv.csv',
    #"structural_fingerprints" : input_dir+prefix+'combination_'+dataset_type+add_in+'_structural_fingerprints.csv',
    #"matrix_fact":'descriptorGeneration/documents-export-2015-12-03/'+dataset_type+'_mat.csv',
    #"neighbours":'descriptorGeneration/documents-export-2015-12-03/'+dataset_type+'_neighbors.csv',
    "gene":final_file_name,
    #"mut":'factors_mutations_50.csv',
    }
    perform_feature_selection = False#True
    grid_search = False # find hyperparameters for the dataset (slow)
    N_JOBS = 3  # the number of jobs to set when running the model or performing grid search
    n_est_grid = 50   # number of estimators when performing grid search
    N_ESTIMATORS = 2000  # manually set the number of estimators in final model
    correlate_scores = False    # Slow - calculates distance between training and test data, and correlates this with difference in prediction of synergy.
    ############################################################################

    X, y, IDs, number_of_structures, stratified_IDs, number_of_occurrences = getXy('train',challenge=2,fileNames)

    if perform_feature_selection == True:
        print "NOTE: FEATURE SELECTION NOT YET IMPLEMENTED IN FINAL MODEL SCRIPT"
        print "Scaling data..."
        # Train a scaler for the features (required for feature selection)
        scaler = preprocessing.MinMaxScaler()
        scaler_fitter = scaler.fit(X)
        # Apply to train data
        X = scaler_fitter.transform(X)
        print "Training an elastic net on the training features..."
        enet = featureRank(X, y)
        features_prev = len(X[0])
        print "\nFeatures before feature selection:",features_prev
        X = enet.transform(X)
        print "Features after feature selection: %d (loss of %d - %.2f %%)"%(len(X[0]), features_prev-len(X[0]), 100.*(float(features_prev)-len(X[0]))/float(features_prev))




    if grid_search == True:
        # Grid search
        print "\nPerforming grid search for hyperparameters. This may take a while..."

        #from scipy.stats import randint as sp_randint
        param_dist = {
                    "max_depth": [None],
                    "min_samples_split": [2,3,4,5],
                    "min_samples_leaf": [3,4,5],
                    "bootstrap": [True],
                    'max_features':['auto', 'sqrt', 'log2', None],
                    "criterion":['gini', 'entropy']
                    }
        estimator = ensemble.RandomForestClassifier(n_estimators =n_est_grid , n_jobs=1, random_state=123)
        #crossval = cross_validation.KFold(stratified_IDs, n_folds=5, shuffle=True, random_state=123)
        crossval = cross_validation.StratifiedKFold(stratified_IDs.ravel(), n_folds=5, shuffle=True, random_state=123)
        param_search = GridSearchCV(estimator, param_grid=param_dist, scoring='f1_weighted', n_jobs=N_JOBS, cv=crossval, verbose=1)
        search_results = param_search.fit(X,y.ravel())
        estimator_params = search_results.best_params_
        print "Forcing parameters to include more trees - N=",N_ESTIMATORS
        estimator_params['n_estimators'] = N_ESTIMATORS
        print "Optimised estimator parameters (save these for future use):",estimator_params
        cPickle.dump(estimator_params,open('optimised_params_classifier.pkl','wb'))
    else:
        print "\nSetting pre-determined parameters for model..."
        #estimator_params = {'max_features': 'auto', 'min_samples_split': 2, 'bootstrap': True, 'max_depth': None, 'min_samples_leaf': 4}    # Used in round1
        #estimator_params = {'max_features': 'auto', 'min_samples_split': 2, 'bootstrap': True, 'max_depth': None, 'min_samples_leaf': 3}   # found after performing feature selection
        #{'bootstrap': False, 'min_samples_leaf': 4, 'min_samples_split': 2, 'criterion': 'gini', 'max_features': 'auto', 'max_depth': None}
        estimator_params = cPickle.load(open('optimised_params_classifier.pkl','rb'))
        print "Estimator parameters:",estimator_params

    cv = cross_validation.StratifiedKFold(stratified_IDs,n_folds=10, shuffle=True, random_state=123)

    for estimator_name in ['RFR']:#['SVR','RFR']:
        #p.figure(figsize=(6,6))

        avg_mcc = 0.
        avg_prec = 0.
        avg_rec = 0.
        avg_f1 = 0.
        avg_acc = 0.
        avg_bal_acc = 0.
        avg_auc = 0.
        avg_spec = 0.
        avg_bal_acc = 0.

        for i, (train_i, test_i) in enumerate(cv):

            test_X, test_y = X[test_i],y[test_i]
            train_X, train_y = X[train_i],y[train_i]

            train_y = train_y.ravel()

            if i == 0:
                print "\nTrain set size"
                print "X data:",np.shape(train_X)
                print "y data:",np.shape(train_y)
                print "Test set size"
                print "X data:",np.shape(test_X)
                print "y data:",np.shape(test_y)

            if estimator_name == 'RFR':
                # Set the parameters
                #print "\nForce-setting number of estimators to",N_ESTIMATORS
                estimator = ensemble.RandomForestClassifier(n_estimators =N_ESTIMATORS, n_jobs=N_JOBS, random_state=123)
                print "Applying best settings..."
                #estimator_params['n_estimators'] = 10
                estimator.set_params(**estimator_params)
            elif estimator_name == 'SVR':
                estimator = svm.SVR(random_state=123)

            print "Training..."
            estimator.fit(train_X, train_y)
            possible_classes = list(estimator.classes_)

            print "Predicting..."
            results_prob = estimator.predict_proba(test_X)
            results_class = estimator.predict(test_X)

            print "Analysing..."
            binariser = preprocessing.LabelBinarizer().fit(y)   # from all possible classes, build a binary matrix based upon labels
            binarised_results = binariser.transform(results_class)
            synergy_binary = binarised_results[:,list(binariser.classes_).index('Synergistic')] # get index from binaiser
            synergy_probabilities = results_prob[:,possible_classes.index('Synergistic')]

            binarised_truths =  binariser.transform(test_y)

            print "\nStats for current split:",i
            print metrics.classification_report(test_y,results_class)
            synergy_mcc = metrics.matthews_corrcoef(binarised_truths[:,list(binariser.classes_).index('Synergistic')],binarised_results[:,list(binariser.classes_).index('Synergistic')])
            synergy_prec = metrics.precision_score(binarised_truths[:,list(binariser.classes_).index('Synergistic')],binarised_results[:,list(binariser.classes_).index('Synergistic')])
            synergy_rec = metrics.recall_score(binarised_truths[:,list(binariser.classes_).index('Synergistic')],binarised_results[:,list(binariser.classes_).index('Synergistic')])
            synergy_f1 = metrics.f1_score(binarised_truths[:,list(binariser.classes_).index('Synergistic')],binarised_results[:,list(binariser.classes_).index('Synergistic')])
            synergy_acc = metrics.accuracy_score(binarised_truths[:,list(binariser.classes_).index('Synergistic')],binarised_results[:,list(binariser.classes_).index('Synergistic')])
            fpr, tpr, thresholds = metrics.roc_curve(binarised_truths[:,list(binariser.classes_).index('Synergistic')],synergy_probabilities)
            synergy_roc_auc = metrics.auc(fpr,tpr)

            cm = metrics.confusion_matrix(binarised_truths[:,list(binariser.classes_).index('Synergistic')],binarised_results[:,list(binariser.classes_).index('Synergistic')])
            TN, FN = float(cm[0][0]),float(cm[1][0])
            TP, FP = float(cm[1][1]), float(cm[0][1])
            specificity = TN/(TN+FN)
            bal_acc = (synergy_prec + specificity) / 2.

            avg_mcc += synergy_mcc
            avg_prec += synergy_prec
            avg_rec += synergy_rec
            avg_f1 += synergy_f1
            avg_acc += synergy_acc
            avg_auc += synergy_roc_auc
            avg_spec += specificity
            avg_bal_acc += bal_acc

            print "\nMCC=%.2f, prec=%.2f, rec=%.2f, F1=%.2f, acc=%.2f, auc=%.2f, bal_acc=%.2f, spec=%.2f"%(synergy_mcc,synergy_prec,synergy_rec,synergy_f1,synergy_acc,synergy_roc_auc,bal_acc,specificity)

        avg_mcc /= len(cv)
        avg_prec /= len(cv)
        avg_rec /= len(cv)
        avg_f1 /= len(cv)
        avg_acc /= len(cv)
        avg_auc /= len(cv)
        avg_spec /= len(cv)
        avg_bal_acc /= len(cv)

        print "\nAverages: MCC=%.2f, prec=%.2f, rec=%.2f, F1=%.2f, acc=%.2f, auc=%.2f, bal_acc=%.2f, spec=%.2f"%(avg_mcc,avg_prec,avg_rec,avg_f1,avg_acc,avg_auc, avg_spec, avg_bal_acc)
        out.write("\n"+final_file_name)
        out.write("\nAverages: MCC=%.2f, prec=%.2f, rec=%.2f, F1=%.2f, acc=%.2f, auc=%.2f, bal_acc=%.2f, spec=%.2f"
                  %(avg_mcc,avg_prec,avg_rec,avg_f1,avg_acc,avg_auc, avg_spec, avg_bal_acc))

print "SD Expressions, models are Done :)"
