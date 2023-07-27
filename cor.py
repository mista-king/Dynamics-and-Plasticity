import os
import time
import pandas as pd
from multiprocessing import Pool, cpu_count


def write_result(batch_res):
    with open('./result_cor.txt', 'a+') as f:
        for res in batch_res:
            f.write("{},{},{}\n".format(res[0], res[1], res[2]))


def error_back(err):
    print(err)


def process_func(i, pt):
    batch_result = []
    print("proceed:{} - task:{}".format(os.getpid(), i))
    for j in range(i+1, len(pt.columns)):
        cor = pt.iloc[:, [i, j]].corr().iat[0, 1]
        batch_result.append((i, j, cor))
    return batch_result


if __name__ == '__main__':
    df = pd.read_csv("/dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4/MT/OvarianCarcinoma_chrMmutation.csv", usecols=[1, 2, 3])
    pt = pd.pivot_table(df, values='frequency', index='mutation', columns='position') 
    pt = pt.fillna(0)
    cols = pt.columns

    cols_df = pd.DataFrame(cols)
    cols_df.to_csv("./cols_index_cor.csv")

    print("CPUcores:{}".format(cpu_count()))
    print("allproceed:{}".format(os.getpid()))

    pool = Pool(cpu_count())
    start = time.time()

    for i in range(len(cols)):
        pool.apply_async(func=process_func, args=(i, pt,), callback=write_result, error_callback=error_back)

    pool.close()
    pool.join()

    end = time.time()
    print("total{:.2f}seconds".format(end-start))
