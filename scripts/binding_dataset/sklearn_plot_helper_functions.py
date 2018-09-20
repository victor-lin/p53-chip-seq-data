import pandas as pd
import numpy as np
import re

from sklearn.svm import SVC
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    matthews_corrcoef,
    roc_auc_score
)
from sklearn.feature_selection import RFE
from xgboost.sklearn import XGBClassifier

import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="white", palette="muted", color_codes=True)


def get_svc_coefficients(df_train, df_test):
    """Return coefficients of SVC model.

    - train SVC with C=0.01 and balanced class weights
    - print classification report
    """
    X_train = df_train.iloc[:, 1:].values
    y_train = df_train.iloc[:, 0].values
    X_test = df_test.iloc[:, 1:].values
    y_test = df_test.iloc[:, 0].values

    clf = SVC(kernel='linear', class_weight='balanced', C=0.01)
    clf.fit(X_train, y_train)

    print classification_report(y_test, clf.predict(X_test))
    print confusion_matrix(y_test, clf.predict(X_test))
    return clf.coef_[0]


def get_svc_rfe_rankings(df_train, df_test):
    """Return feature rankings of SVC model determined by recursive feature elimination.

    - train SVC with C=0.01 and balanced class weights
    - print classification report
    """
    X_train = df_train.iloc[:, 1:].values
    y_train = df_train.iloc[:, 0].values
    X_test = df_test.iloc[:, 1:].values
    y_test = df_test.iloc[:, 0].values

    estimator = SVC(kernel='linear', class_weight='balanced', C=0.01)
    selector = RFE(estimator, n_features_to_select=1, step=1)
    selector = selector.fit(X_train, y_train)

    print classification_report(y_test, selector.predict(X_test))
    print confusion_matrix(y_test, selector.predict(X_test))
    return selector.ranking_


def get_xgboost_importances(df_train, df_test):
    """Return feature importances of trained XBGC model.

    - train XGBC with 20:1 negative to positive ratio
    - print classification report
    """
    X_train = df_train.iloc[:, 1:].values
    y_train = df_train.iloc[:, 0].values
    X_test = df_test.iloc[:, 1:].values
    y_test = df_test.iloc[:, 0].values

    clf = XGBClassifier(nthreads=-1, scale_pos_weight=20)
    clf.fit(X_train, y_train)

    print classification_report(y_test, clf.predict(X_test))
    print confusion_matrix(y_test, clf.predict(X_test))
    return clf.feature_importances_


def get_abs_feature_weights_dict(weights, names):
    """Return dict of {feature name: |weight|}"""
    weights = [abs(x) for x in weights]
    return dict(zip(names, weights))


def get_motif_weights_plot(weights, names, xmax, title=None):
    """Return plot of motif features with absolute values of weights."""
    feature_weights = get_abs_feature_weights_dict(weights, names)

    motif_weights = {x: feature_weights[x]
                     for x in feature_weights
                     if 'P53match' in x}

    df_motif_ranks = pd.DataFrame.from_dict(motif_weights, orient='index')
    df_motif_ranks['motif'] = [re.sub('P53match_(score|count)(_(sum|max))*_', '', f)
                               for f in df_motif_ranks.index]
    df_motif_ranks['type'] = [re.match(r'P53match_(score_max|score_sum|count).*', f).groups()[0]
                              for f in df_motif_ranks.index]
    df_motif_ranks.reset_index(drop=True, inplace=True)

    all_motifs = list(set(df_motif_ranks['motif']))
    all_types = list(set(df_motif_ranks['type']))

    # fig, ax = plt.subplots(figsize=(10, 5))
    fig, ax = plt.subplots()
    ind = np.arange(len(set(df_motif_ranks['motif'])))
    width = 0.25
    colors = ['r', 'y', 'g']

    type_bars = []
    for i, t in enumerate(all_types):
        data = [df_motif_ranks[(df_motif_ranks['type'] == t) & (df_motif_ranks['motif'] == motif)][0].values[0] for motif in all_motifs]
        p = ax.barh(ind - width + width * i, data, width, color=colors[i])
        type_bars.append(p)

    ax.set_yticks(ind)
    ax.set_yticklabels(all_motifs)
    ax.set_xlim(0, xmax)
    ax.set_xlabel('weight')
    ax.legend([x[0] for x in type_bars], all_types)

    # ax.autoscale_view()
    if title:
        plt.title(title)
    plt.grid()
    return plt


def get_all_weights_plot(weights, names, xmax, title=None):
    """Return plot of top 20 important features determined by absolute value of weight."""
    imp = [abs(x) for x in weights]
    imp, names = zip(*sorted(zip(imp, names), reverse=True)[:20][::-1])
    plt.barh(range(len(names)), imp, align='center')
    plt.yticks(range(len(names)), names)
    if title:
        plt.title(title)
    plt.xlim(0, xmax)
    plt.xlabel('weight')
    return plt


def get_gridsearch_scores_plot(grid_clf, C_values):
    """Return plot of validation scores from a gridsearch over values for C parameter of SVC."""
    scores = [x.mean_validation_score for x in grid_clf.grid_scores_]
    # score_stds = [np.std(x.cv_validation_scores) for x in grid_clf.grid_scores_]

    plt.plot(C_values, scores)
    # plt.legend()
    plt.xlabel('C')
    plt.xticks(C_values)
    plt.xscale('log', basex=2)
    plt.ylabel('Mean score')
    # plt.ylim(0.3,0.7)
    return plt


def get_svc_param_plot(df_train, df_test):
    """Return plot of MCC and ROCAUC scoring parameters against values for C parameter of SVC."""
    X_train = df_train.iloc[:, 1:].values
    y_train = df_train.iloc[:, 0].values
    X_test = df_test.iloc[:, 1:].values
    y_test = df_test.iloc[:, 0].values

    C_values = np.logspace(-9, -6, num=4, base=2)
    mccs = []
    roc_aucs = []

    for C in C_values:
        clf = SVC(kernel='linear', class_weight='balanced', C=C)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        names = df_train.columns[1:]
        imp = [abs(x) for x in clf.coef_[0]]
        imp, names = zip(*sorted(zip(imp, names), reverse=True)[:20])

        mcc = matthews_corrcoef(y_test, y_pred)
        roc_auc = roc_auc_score(y_test, y_pred)
        mccs.append(mcc)
        roc_aucs.append(roc_auc)

    plt.plot(C_values, mccs)
    plt.plot(C_values, roc_aucs)
    plt.xscale('log', basex=2)
    return plt
