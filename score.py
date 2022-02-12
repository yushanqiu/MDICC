import pandas as pd
import numpy as np
from sklearn.cluster import KMeans

def MDICCscore(x,k,path,label):
    k = int(k)
    kmeans = KMeans(n_clusters=k,random_state = 0).fit(x)
    py_result = kmeans.predict(x)

    true_label = pd.read_csv(path)
    true_label = true_label[label].tolist()
    true_label1 = np.array(true_label)
    a = np.array(np.where(np.isnan(true_label1)))
    r, c = a.shape
    if c > 0:
        true_label1 = np.delete(true_label1, a)
        py_result = np.delete(py_result, a)
    true_label = true_label1.tolist()
    test_label = py_result.tolist()

    from sklearn.metrics.cluster import rand_score
    RI = rand_score(true_label, test_label)

    from sklearn.metrics import adjusted_rand_score
    ARI = adjusted_rand_score(true_label, test_label)

    from sklearn.metrics.cluster import normalized_mutual_info_score
    NMI = normalized_mutual_info_score(true_label, test_label)

    from sklearn.metrics import accuracy_score
    accu = accuracy_score(np.array(true_label), np.array(test_label))

    from sklearn.metrics import f1_score
    F1 = f1_score(np.array(true_label), np.array(test_label), average='micro')


    level = np.array((RI,ARI,NMI,accu,F1), dtype = "float16")

    return level