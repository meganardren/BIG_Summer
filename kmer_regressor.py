#Code adapted from Mudi Yang
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.linear_model

def kmer_regressor(tile_num):
    #kmer_size used to name trained model and to name processed data
    kmer_size = '6mer'
    Data_path = "/u/home/m/mardren/scratch/KmerML"
    #read in data
    label_tile = pd.read_pickle('%s/%s_sequences_label.pkl'%(Data_path, kmer_size))
    encoding_tile = pd.read_pickle('%s/%s_sequences_encoding.pkl'%(Data_path, kmer_size))
#     #subset by selected tile
#     label_tile = sequence_label[[tile_num]]
#     encoding_tile = sequence_encoding[[tile_num]]
    #reformat data
    label_tile =label_tile.unstack().reset_index()
    label_tile = label_tile.rename(columns={0:'label'})
    encoding_tile =encoding_tile.unstack().reset_index()
    encoding_tile = encoding_tile.rename(columns={0:'encoding'})
    encoding_tile = encoding_tile.set_index('region_id')
    label_tile = label_tile.set_index('region_id')
    encoding_tile = encoding_tile.drop(columns=['tile_num'])
    all_data = pd.concat([label_tile,encoding_tile], axis=1)
    all_data = all_data.reset_index()
    all_data['chrom'] = all_data['region_id'].str.split('_').str[3]

    #setting aside chr1 and chr2 for testing data
    train_chroms = [f'chr{i}' for i in range(1,23,2)]
    test_chroms = [f'chr{i}' for i in range(2,23,2)]
    chrX = all_data[all_data["chrom"]==('chrX')]
    test_set = all_data[all_data["chrom"].isin(test_chroms)]
    train_set = all_data[all_data["chrom"].isin(train_chroms)]
    train_set = train_set.append(chrX)

    test_set = test_set.drop(columns=['chrom'])
    train_set = train_set.drop(columns=['chrom'])
    all_data = all_data.drop(columns=['chrom'])
    

    #format the training data
    labels = train_set['label'].values
    encodings = train_set['encoding'].values.tolist()
    encodings = np.asarray(encodings)
    labels = np.reshape(labels,(258633 ,))
    #format the testing data
    test_labels = test_set['label'].values
    test_encodings = test_set['encoding'].values.tolist()
    test_encodings = np.asarray(test_encodings)
    test_labels = np.reshape(test_labels,(228687,))

    encodings_df = pd.DataFrame(encodings)
    encodings_df.fillna(encodings_df.mean(), inplace=True)
    encodings = encodings_df.to_numpy()
    labels_df = pd.DataFrame(labels)
    labels_df.fillna(labels_df.mean(), inplace=True)
    labels = labels_df.to_numpy()
    labels = np.reshape(labels,(258633 ,))

    test_encodings_df = pd.DataFrame(test_encodings)
    test_encodings_df.fillna(test_encodings_df.mean(), inplace=True)
    test_encodings = test_encodings_df.to_numpy()
    test_labels_df = pd.DataFrame(test_labels)
    test_labels_df.fillna(test_labels_df.mean(), inplace=True)
    test_labels = test_labels_df.to_numpy()
    test_labels = np.reshape(test_labels,(228687,))

    #format all data
    data_labels = all_data['label'].values
    data_encodings = all_data['encoding'].values.tolist()
    data_encodings = np.asarray(data_encodings)
    data_labels = np.reshape(data_labels,(487320,))

    data_encodings_df = pd.DataFrame(data_encodings)
    data_encodings_df.fillna(data_encodings_df.mean(), inplace=True)
    data_encodings = data_encodings_df.to_numpy()
    data_labels_df = pd.DataFrame(data_labels)
    data_labels_df.fillna(data_labels_df.mean(), inplace=True)
    data_labels = data_labels_df.to_numpy()

    data_labels = np.reshape(data_labels,(487320,))

    #train a linear regressor on the training dataset. 
    from sklearn.linear_model import LinearRegression

    k_mer_regressor = LinearRegression()
    k_mer_regressor.fit(encodings, labels)

    from sklearn.metrics import mean_squared_error
    test_predictions = k_mer_regressor.predict(test_encodings)
    test_mse = mean_squared_error(test_predictions, test_labels)
    test_rmse = np.sqrt(test_mse)

    #display
    caption = "correlation pearson r= %s"%(pd.Series(test_predictions).corr(pd.Series(test_labels)))
    mse = "root mean squared error = %s"%test_rmse
    plt.scatter(
        x=test_predictions,
        y=test_labels,
        edgecolors='w'
    )

    plt.title("%s Regressor on Even testing data"%kmer_size)
    plt.xlabel("test_predictions \n\n %s"%caption)
    plt.ylabel("test_labels")

    plt.savefig('6mer_eventesting')

    #get predictions for the whole tile
    predictions = k_mer_regressor.predict(data_encodings)
    all_data['prediction'] = predictions
    #save results to use in sequence model
    all_data = all_data.drop(columns=['encoding'])
    #split into kmer_prediction pkl and sharpr_score pkl
    kmer_prediction = all_data.drop(columns=['label'])
    kmer_prediction = kmer_prediction.pivot(index='region_id', columns='tile_num', values='prediction')
    sharpr_score = all_data.drop(columns=['prediction'])
    sharpr_score = sharpr_score.pivot(index='region_id', columns='tile_num', values='label')
    tile_num = tile_num+1
    kmer_prediction.to_pickle('%s/%s_prediction_tile_evenchrtesting.pkl'%(Data_path, kmer_size))
    sharpr_score.to_pickle('%s/%s_sharpr_score_tile_evenchrtesting.pkl'%(Data_path, kmer_size))

if __name__ == "__main__":
    tile_num = 15
    kmer_regressor(tile_num)
