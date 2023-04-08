"""
## Classifiers

Non supervised classifiers (random forest, k-nearest neighbors, neural
networks) for predicting the methylation class.
"""

# start_external_modules
import logging
import time
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
# end_external_modules

# start_internal_modules
from nanodip.config import (
    ENDING,
    NANODIP_REPORTS,
)
from nanodip.data import (
    Reference,
    Sample,
    get_sample_methylation,
    reference_methylation_from_index,
)
from nanodip.utils import (
    composite_path,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)


def evaluate_clf(clf, x_sample, X_test, y_test):
    """Calculates classifier accuracy and predicts sample value by
    returning the most probable values as formatted string.

    Args:
        x_sample: Feature list of sample.
        X_test: Features of test data.
        y_test: Classes of test data.

    Returns:
        Formatted string containing most probable classes of x_sample
        and classifier accuracy.
    """
    y_predict = clf.predict(X_test)
    # Fraction of correctly classified test samples.
    accuracy = accuracy_score(y_test, y_predict)
    prob = clf.predict_proba([x_sample])[0]
    prob_per_class = list(zip(prob, clf.classes_))
    prob_per_class.sort(reverse=True)
    result = (
        "Evaluation of %s\n"
        "Classifier accuracy: %.2f %%\n"
        "Classifier probability per class:\n"
    ) % (clf, 100*accuracy)
    for p, c in prob_per_class[:10]:
        result += "%-16s : %5.2f %%\n" % (c, 100*p)
    return result

def training_test_data(sample, reference):
    """Takes the reference data that overlaps with the sample CpGs and
    splits it into training data and test data.

    Args:
        sample: Sample to determine CpG overlap.
        reference: Comparison reference.

    Returns:
        X_train, X_test, y_train, y_test: Split training and test
            data pairs.
    """
    X = reference_methylation_from_index(
        reference.specimens_index, sample.cpg_overlap_index
    )
    y = reference.methylation_class
    return train_test_split(X, y, test_size=0.2)

def fit_and_evaluate_classifiers(sample_name, reference_name):
    """Uses non supervised machine learning classifiers (random forest,
    k-nearest neighbors, neural networks) to predict the methylation
    class of the sample. Output will be written to disk.

    Args:
        sample_name: Name of sample to analyze.
        reference_name: Name of reference that is used to train classifiers.
    """
    sample = Sample(sample_name)
    reference = Reference(reference_name)
    # Define training/test/sample data.
    X_train, X_test, y_train, y_test = training_test_data(sample, reference)
    x_sample = get_sample_methylation(sample)
    # Define classifier models.
    rf_clf = RandomForestClassifier(
        n_estimators=150,
        n_jobs=-1,
    )
    knn_clf = KNeighborsClassifier(
        n_neighbors=5,
        weights="distance",
    )
    nn_clf = MLPClassifier(
        verbose=True,
    )
    # SVM are very time consuming.
    # from sklearn.svm import SVC
    # svm_clf = SVC(
        # kernel="linear",
        # probability=True,
        # verbose=True,
    # )
    clfs = [rf_clf, knn_clf, nn_clf]
    output_file = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["clf_txt"],
    )
    # Clean file.
    with open(output_file, "w") as f:
        f.write("")
    # Stop if there is no data to fit.
    if len(x_sample) == 0:
        with open(output_file, "a") as f:
            f.write(f"No data to fit.\n")
        return
    # Otherwise train classifiers and evaluate.
    for clf in clfs:
        with open(output_file, "a") as f:
            f.write(f"Start training {clf}.\n")
        start = time.time()
        clf.fit(X_train, y_train)
        evaluation = evaluate_clf(clf, x_sample, X_test, y_test)
        passed_time = time.time() - start
        with open(output_file, "a") as f:
            f.write("Time used for classification: %.2f s\n" % passed_time)
            f.write(evaluation + "\n")
