import pandas as pd
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

folderName = 'input_folder/'

### Generating and training 5 random forest models ###
for i in range(1,6):
    modelName = folderName + 'model_' + str(i) +'.json'

    #loading the data for the current training group
    df = pd.read_table(folderName + 'train_' + str(i) + '.txt')

    #setting the  data and the labels for the training group
    y_train = df.is_true_somatic
    X_train = df
    X_train = X_train.iloc[:, :-1]

    #generating and training the random forest model on the training group
    clf = RandomForestClassifier(n_estimators=50, random_state=42)
    clf.fit(X_train,y_train)
    pickle.dump(clf,open(modelName,'wb'))

    # loading the data for the current validation group
    df_val = pd.read_table(folderName + 'val_' + str(i) + '.txt')

    # setting the  data and the labels for the validation group
    y_val = df_val.is_true_somatic
    X_val = df_val
    X_val = X_val.iloc[:, :-1]

    #predicting calssificatoin for validation group
    y_pred=clf.predict(X_val)

    # estimating models' performance on validation set
    print(metrics.confusion_matrix(y_val, y_pred))
    print('first precision ', precision_score(y_val, y_pred))
    print('first recall ' , recall_score(y_val, y_pred))


### testing models perforance on test set ###
test = pd.read_table(folderName + 'test.txt')
test_data = test.iloc[:, :-1]
real = test.is_true_somatic

#for majority voting
sum_pred = np.zeros(len(real))

# using all 5 trained models for prediction #
for i in range(1,6):

    #loading the trained model
    modelName1 = folderName + 'model_' + str(i) + '.json'
    model = pickle.load(open(modelName1,'rb'))

    #making the predictions on the test group
    pred1 = model.predict(test_data)
    sum_pred += pred1

#actual prediction is based on majority voting
i=np.where(sum_pred>2)
i=pd.DataFrame(i)
i=i.T
final_pred = np.zeros(len(sum_pred))
np.put(final_pred,i,1)

# All predicted somatic mutations are tested via RNA-MuTect filtering steps.
# Mutations that are filtered out will not be considered as somatic mutations.

# estimating the performance on the test set after RNA-MuTect filtering steps
print('final precision ', precision_score(real, final_pred))
print('final recall ' , recall_score(real, final_pred))