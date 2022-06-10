"""
## Classifiers

Non supervised classifiers (Random forest, k-nearest neighbors, neural
networks, support vector machines) for predicting the methylation class.
"""

import time

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC

from config import (
    ENDING,
    NANODIP_REPORTS,
)
from data import (
    Reference,
    Sample,
    get_sample_methylation,
    reference_methylation_from_index,
)
from utils import (
    composite_path,
)

def evaluate_clf(clf, x_sample, X_test, y_test):
    """Calculates classifier accuracy and predicts sample value by
    returning the most probable values as formatted string.
    """
    y_predict = clf.predict(X_test)
    # Fraction of correctly classified test samples.
    accuracy = accuracy_score(y_test, y_predict)
    prob = clf.predict_proba([x_sample])[0]
    prob_per_class = [(p, c) for p, c in zip(prob, clf.classes_)]
    prob_per_class.sort(reverse=True)
    result = (
        "Evaluation of %s\n"
        "Classifier accuracy: %.2f %%\n"
        "Classifier probability per class:\n"
    ) % (clf, 100*accuracy)
    for i in range(10):
        result += (
            "%-16s : %5.2f %%\n" % (
                prob_per_class[i][1],
                100*prob_per_class[i][0],
            )
        )
    return result

def training_test_data(sample, reference):
    """Takes the reference data that overlaps with the sample CpGs and
    splits it into training data and validation data.
    """
    sample.set_cpg_overlap(reference)
    X = reference_methylation_from_index(
        reference.specimens_index, sample.cpg_overlap_index
    )
    y = reference.methylation_class
    return train_test_split(X, y, test_size=0.2)

def fit_and_evaluate_classifiers(sample_name, reference_name):
    """Uses non supervised machine learning classifiers (Random forest,
    k-nearest neighbors, neural networks, support vector machines)
    to evaluate methylation class of sample.

    Args:
        sample_name: Name of sample to analyze.
        reference: Name of reference that is used to train classifiers.
    Returns:
        TODO
    """
    sample = Sample(sample_name)
    reference = Reference(reference_name)
    # Define training/test/sample data.
    X_train, X_test, y_train, y_test = training_test_data(sample, reference)
    x_sample = get_sample_methylation(sample, reference)
    # Define classifier models.
    rf_clf = RandomForestClassifier(
        n_estimators=150,
        n_jobs=-1,
        random_state=1234,
    )
    knn_clf = KNeighborsClassifier(
        n_neighbors=5,
        weights="distance",
    )
    nn_clf = MLPClassifier()
    svm_linear_clf = SVC(
        kernel="linear",
        probability=True,
    )
    clfs = [rf_clf, knn_clf, nn_clf, svm_linear_clf]
    # clfs = [rf_clf, knn_clf, nn_clf]
    output_file = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["clf"],
    )
    with open(output_file, "w") as f:
        f.write("")
    # Train classifiers and evaluate.
    for clf in clfs:
        with open(output_file, "a") as f:
            f.write(f"Start training {clf}.\n")
        start = time.time()
        clf.fit(X_train, y_train)
        evaluation = evaluate_clf(clf, x_sample, X_test, y_test)
        passed_time = time.time() - start
        with open(output_file, "a") as f:
            f.write(f"Time used for classification: %.2f s\n" % passed_time)
            f.write(evaluation + "\n")
