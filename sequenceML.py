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

    #check correlation
    sharpr_15_array = np.asarray(sharpr_15)
    sharpr_15_array = sharpr_15_array.flatten()
    conv_15_array = np.asarray(conv_15)
    conv_15_array = conv_15_array.flatten()
    deepfact_15_array = np.asarray(deepfact_15)
    deepfact_15_array = deepfact_15_array.flatten()
    kmer_15_array = np.asarray(kmer_15)
    kmer_15_array = kmer_15_array.flatten()
    gkm_15_array = np.asarray(gkm_15)
    gkm_15_array = gkm_15_array.flatten()
    print('corr of SHARPR score and dragoNN ConvModel:')
    print(pd.Series(sharpr_15_array).corr(pd.Series(conv_15_array)))
    print('corr of SHARPR score and gkmSVMModel:')
    print(pd.Series(sharpr_15_array).corr(pd.Series(gkm_15_array)))
    print('corr of SHARPR score and dragoNN DeepFactorizedModel:')
    print(pd.Series(sharpr_15_array).corr(pd.Series(deepfact_15_array)))
    print('corr of SHARPR score and 6merModel:')
    print(pd.Series(sharpr_15_array).corr(pd.Series(kmer_15_array)))
    

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
    data = data.drop(columns=['chrom'])

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


    #format data
    data_labels = data['sharpr_%s'%tile_num].values
    data = data.drop(columns=['region_id','sharpr_%s'%tile_num])
    data_encodings = data.values.tolist()
    data_encodings = np.asarray(data_encodings)
    data_encodings_df = pd.DataFrame(data_encodings)
    data_encodings_df.fillna(data_encodings_df.mean(), inplace=True)
    data_encodings = data_encodings_df.to_numpy()
    data_labels = np.reshape(data_labels,(15720,))

    #train a linear regressor on the training dataset. 
    from sklearn.linear_model import LinearRegression

    sequence_regressor = LinearRegression()
    sequence_regressor.fit(encodings, labels)

    #compute mse
    #     from sklearn.metrics import mean_squared_error
    #     test_predictions = sequence_regressor.predict(test_encodings)
    #     test_mse = mean_squared_error(predictions, data_labels)
    #     test_rmse = np.sqrt(test_mse)
    from sklearn.metrics import mean_squared_error
    predictions = sequence_regressor.predict(data_encodings)
    mse = mean_squared_error(predictions, data_labels)
    rmse = np.sqrt(mse)

    data_encodings=data_encodings.flatten()
    data__labels=data_labels.flatten()
    #visualize the results
    corr = pd.Series(predictions).corr(pd.Series(data_labels))

    caption = "correlation pearson r= %s"%(corr)
    rmse = "root mean squared error = %s"%rmse


    plt.scatter(
        x=predictions,
        y=data_labels,
        edgecolors='w'
    )
    tile_num = int(tile_num)
    tile_num = tile_num + 1

    plt.title("Tile %s Sequence Regressor"%(tile_num))
    plt.xlabel("test_predictions \n\n %s \n\n %s"%(caption,rmse))
    plt.ylabel("test_labels")

    plt.savefig("sequence_predictions_tile%s"%tile_num)

    print("Regressor coefficients(ConvModel,gkmSVMModel,DeepFactorizedModel,6merModel):")
    print(sequence_regressor.coef_)

#     fn = open("sequence_predictions_tile%s.txt"%tile_num,"w")
#     np.savetxt(fn,predictions)
#     fn.close()


if __name__ == "__main__":
    # take tile_num in as a parameter
    tile_num = '15'
    sequence_file(tile_num)
