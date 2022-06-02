import pandas as pd
import numpy as np
from config import *
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
import pickle
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
import os
import matplotlib.pyplot as plt


def prepare_train_test():
    os.system(f"rm -f {cancer_dir}/*_set.csv")
    final_cols_to_use = features + ['Chromosome', 'position', 'Tumor_Sample_Barcode', 'judgement', 'is_real_keep',
                                    'ref_allele', 'alt_allele']
    full_data = pd.read_csv(os.path.abspath(cancer_dir + "/all_vars_after_DNA_pon.csv"),
                            usecols=final_cols_to_use,dtype={'Chromosome': 'str'})
    full_data = full_data.drop_duplicates()

    # making sure all features are numeric and not nan
    for feature in features:
        full_data[feature] = pd.to_numeric(full_data[feature])
        cur_mean = np.nanmean(full_data[feature], axis=0)
        full_data[feature] = np.where(full_data[feature].isna(), cur_mean, full_data[feature])
    # split to train-test
    sample_names = full_data["Tumor_Sample_Barcode"].drop_duplicates()
    train_sample_names = set(np.random.choice(sample_names, num_train_samples, replace=False))
    test_sample_names = set(sample_names) - train_sample_names
    assert (not any(name in train_sample_names for name in test_sample_names))
    assert (not any(name in test_sample_names for name in train_sample_names))
    test_set = full_data[full_data["Tumor_Sample_Barcode"].isin(
        test_sample_names)]
    assert (test_set['Tumor_Sample_Barcode'].drop_duplicates().shape[0] == (
            full_data['Tumor_Sample_Barcode'].drop_duplicates().shape[0] - num_train_samples))
    train_set = full_data[full_data["Tumor_Sample_Barcode"].isin(train_sample_names)]
    assert (train_set['Tumor_Sample_Barcode'].drop_duplicates().shape[0] == num_train_samples)
    # removing noise variants, keeping only germline (REJECT) or somatic (KEEP and is_real_keep=1)
    train_set = train_set[~(train_set['judgement'] == 'KEEP') | ~(
            train_set['is_real_keep'] == 0)]
    return train_set, test_set


def train_models(train_set):
    if not os.path.isdir(models_path):
        os.mkdir(models_path)
    else:
        os.system(f"rm -f {models_path}/*")  # clear old models
    if not os.path.isdir(results_path):
        os.mkdir(results_path)
    else:
        os.system(f"rm -f {results_path}/train*")  # clear old results
    cv = KFold(n_splits=num_folds)
    train_samples = train_set['Tumor_Sample_Barcode'].drop_duplicates().tolist()
    precision_scores, recall_scores, validation_samples, folds = [], [], [], []
    with open(os.path.abspath(results_path + "/train_results.txt"), 'w') as results:
        results.write("Training scores\n\n")
        # cross validation of the train samples
        for (train, validation), i in zip(cv.split(train_samples), range(1, num_folds + 1)):
            folds += [i] * len(validation)
            cur_train_samples = np.array(train_samples)[train]
            cur_validation_samples = np.array(train_samples)[validation]
            validation_samples += cur_validation_samples.tolist()
            assert (not any(name in cur_train_samples for name in cur_validation_samples))
            assert (not any(name in cur_validation_samples for name in cur_train_samples))
            cur_train = train_set.loc[train_set["Tumor_Sample_Barcode"].isin(cur_train_samples), :]
            cur_validation = train_set.loc[train_set["Tumor_Sample_Barcode"].isin(cur_validation_samples), :].copy(
                deep=True)

            # train current model
            clf = RandomForestClassifier(n_estimators=50, random_state=42)
            clf.fit(cur_train[features], cur_train['is_real_keep'])
            model_path = os.path.abspath(f"{models_path}model_{str(i)}.json")
            pickle.dump(clf, open(model_path, 'wb'))

            # validation
            y_pred = clf.predict(cur_validation[features])
            y_true = cur_validation['is_real_keep']

            # calc model scores
            precision = precision_score(y_true, y_pred)
            recall = recall_score(y_true, y_pred)
            results.write(f"Model {i}:\n")
            results.write(f'\tPrecision score: {precision}\n')
            results.write(f'\tRecall score: {recall}\n\n')

            # calc scores per sample
            cur_validation.loc[:, 'pred'] = y_pred
            for sample in cur_validation_samples:
                cur_pred = cur_validation[cur_validation["Tumor_Sample_Barcode"] == sample]['pred']
                cur_truth = cur_validation[cur_validation["Tumor_Sample_Barcode"] == sample]['is_real_keep']
                precision_scores.append(precision_score(cur_truth, cur_pred))
                recall_scores.append(recall_score(cur_truth, cur_pred))
            del cur_truth
            del cur_pred
            del cur_train
            del cur_validation
            del cur_train_samples
            del cur_validation_samples
    # save scores per sample file
    scores = pd.DataFrame(
        {'Fold': folds, 'Sample': validation_samples, 'Precision': precision_scores, 'Recall': recall_scores})
    scores.to_csv(os.path.abspath(results_path + "/train_scores_per_sample.csv"), index=False)
    assert (scores.shape[0] == scores['Sample'].drop_duplicates().shape[0])
    # save results boxplot
    fig = plt.figure()
    fig.suptitle(f'Validation Set (n={num_train_samples})', fontsize=20)
    scores[['Precision', 'Recall']].boxplot(grid=False)
    fig.savefig(os.path.abspath(results_path + "/train_boxplot"))
