import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.linear_model
from sklearn.linear_model import LinearRegression
import joblib
import sys

def sequence_file(tile_num):
    #create datapaths 
    Data_path = "/u/home/m/mardren/scratch/SequenceML"
    gkm_path = "%s/features.gkSVM.HepG2.tsv.gz"%Data_path
    conv_path       = "%s/features.dragoNN_ConvModel.HepG2.SV40P.Rep1.tsv.gz"%Data_path
    deepfact_path       = "%s/features.dragoNN_DeepFactorizedModel.HepG2.minP.Rep1.tsv.gz"%Data_path
    kmer_path = "%s/6mer_prediction.pkl"%Data_path
    sharpr_path = "%s/6mer_label.pkl"%Data_path


    #open data from csv into dataframes
    gkm = pd.read_csv(gkm_path, header = 0, index_col = 0, sep = '\t')
    conv  = pd.read_csv(conv_path, header = 0, index_col = 0, sep = '\t')
    deepfact = pd.read_csv(deepfact_path, header = 0, index_col = 0, sep = '\t')
    kmer = pd.read_pickle('%s'%(kmer_path))
    #subset by selected tile
    sharpr = pd.read_pickle('%s'%(sharpr_path))
    gkm_15 = gkm[['feat_gksvm_%s'%tile_num]]
    gkm_15 = gkm_15.rename(columns = {'feat_gksvm_%s'%tile_num:'gkm_%s'%tile_num})
    conv_15 = conv[[tile_num]]
    conv_15 = conv_15.rename(columns = {tile_num:'conv_%s'%tile_num})
    deepfact_15 = deepfact[[tile_num]]
    deepfact_15 = deepfact_15.rename(columns = {tile_num:'deepfact_%s'%tile_num})
    tile_num = int(tile_num)
    kmer_15 = kmer[[tile_num]]
    kmer_15 = kmer_15.rename(columns = {tile_num:'kmer_%s'%tile_num})
    sharpr_15 = sharpr[[tile_num]]
    sharpr_15 = sharpr_15.rename(columns = {tile_num:'sharpr_%s'%tile_num})

    #create data dataframe
    data = pd.concat([conv_15, gkm_15, deepfact_15, kmer_15, sharpr_15], axis=1)

    #reset index to be able to access region_id
    data = data.reset_index()

    #split region_id to create chrom column
    data['chrom'] = data['region_id'].str.split('_').str[3]

    #setting aside chr21 for testing data
    test_set = data[data["chrom"] == 'chr21']
    #removing chr21 from training data
    train_set = data[data["chrom"] != 'chr21']
    #dropping the chrom column
    test_set = test_set.drop(columns=['chrom'])
    train_set = train_set.drop(columns=['chrom'])
    
    #format training data
    labels = train_set['sharpr_%s'%tile_num].values
    train_set = train_set.drop(columns=['region_id','sharpr_%s'%tile_num])
    encodings = train_set.values.tolist()
    encodings = np.asarray(encodings)
    encodings_df = pd.DataFrame(encodings)
    encodings_df.fillna(encodings_df.mean(), inplace=True)
    encodings = encodings_df.to_numpy()
    labels = np.reshape(labels,(15537,))
    
    #format the testing data
    test_labels = test_set['sharpr_%s'%tile_num].values
    test_set = test_set.drop(columns=['region_id','sharpr_%s'%tile_num])
    test_encodings = test_set.values.tolist()
    test_encodings = np.asarray(test_encodings)
    test_encodings_df = pd.DataFrame(test_encodings)
    test_encodings_df.fillna(test_encodings_df.mean(), inplace=True)
    test_encodings = test_encodings_df.to_numpy()
    test_labels = np.reshape(test_labels,(183,))
    
    #train a linear regressor on the training dataset. 
    from sklearn.linear_model import LinearRegression

    sequence_regressor = LinearRegression()
    sequence_regressor.fit(encodings, labels)
    
    #compute mse
    from sklearn.metrics import mean_squared_error
    test_predictions = sequence_regressor.predict(test_encodings)
    test_mse = mean_squared_error(test_predictions, test_labels)
    test_rmse = np.sqrt(test_mse)
    
    test_encodings=test_encodings.flatten()
    test__labels=test_labels.flatten()
    #visualize the results
    corr = pd.Series(test_predictions).corr(pd.Series(test_labels))
    caption = "correlation pearson r= %s"%(corr)
    rmse = "root mean squared error = %s"%test_rmse


#     plt.scatter(
#         x=test_predictions,
#         y=test_labels,
#         edgecolors='w'
#     )
#     tile_num = int(tile_num)
#     tile_num = tile_num + 1

#     plt.title("Tile %s chr 21 Sequence Regressor"%(tile_num))
#     plt.xlabel("test_predictions \n\n %s \n\n %s"%(caption,rmse))
#     plt.ylabel("test_labels")

#     plt.savefig("sequence_predictions_tile%s"%tile_num)
    
#     fn = open("sequence_predictions_tile%s.txt"%tile_num,"w")
#     np.savetxt(fn,test_predictions)
#     fn.close()

    return corr

if __name__ == "__main__":
    # take tile_num in as a parameter
    #tile_num = sys.argv[1]
    corr = []
    for tile_num in range(30):
        tile_num = str(tile_num)
        corr_tile = sequence_file(tile_num)
        corr.append(corr_tile)
    
    
    tiles = list(range(1, 31))
     plt.plot(
        tiles,
        corr
    )

    plt.title("Correlation by Tile")
    plt.xlabel("tile")
    plt.ylabel("correlation")

    plt.show("sequence_predictions_tile")
