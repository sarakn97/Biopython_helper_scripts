# Classification Tree in scikit-learn
# Import DecisionTreeClassifier
from sklearn.tree import DecisionTreeClassifier
# Import train_test_split
from sklearn.model_selection import train_test_split
# Import accuracy_score
from sklearn.metrics import accuracy_score

# Split dataset into 80% train, 20% test
# set straify to y in order for the train and test sets to have the same proportion of class labels as the unsplit dataset.
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, stratify = y, random_state = 1)

# Instantiate dt - random state is set to 1 for reproducibility
dt = DecisionTreeClassifier(max_depth= 2, random_state= 1)

#Fit dt to the training set
dt.fit(X_train, y_train)

# Predict the test set labels
y_pred = dt.predict(X_test)
# Evaluate the test-set accuracy
accuracy_score(y_test, y_pred)

# Dec