from sklearn.cluster import KMeans

def MDICClabel(x,k):
    k = int(k)
    kmeans = KMeans(n_clusters=k,random_state = 0).fit(x)
    py_result = kmeans.predict(x)
    test_label = py_result.tolist()
    return test_label