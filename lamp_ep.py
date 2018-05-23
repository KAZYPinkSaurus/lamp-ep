import sys
import math
import pandas as pd
import numpy as np

class Lamp_ep():
    def __init__(self,aA,aQ,aMax_min_sup,aDataset):
        """
        a in (0,1)

        Parameters
        -------
        aMax_min_sup : int
            min_supの最大値

        aDataset : DataFrame
            カラム名指定(['pattern','Nep','Ne'])

        """
        if(aA >= 1 or aA <= 0):
            print("Error:a=!(0~1)")
            exit()
        self.mA = aA
        self.mMax_min_sup = aMax_min_sup
        self.mQ = aQ
        self.mDataset = aDataset

    def __get_min_tau(self):
        """
        τの一番小さいやつを得る
        """
        for i in range(1,self.mMax_min_sup):
            print("i="+str(i)+" "+str(self.mA**i)+"<="+str(self.mQ/len(self.mining_eps_alg(i).index)))
            if self.mA**i <= self.mQ/len(self.mining_eps_alg(i).index):
                return i
        print("τが得られなかった")
        exit()


    def __get_corrected_pvs(self, aTau):
        """
        選ばれたτでパターン抽出を行い全てで多重検定をい,パターンとp値を返す
        """
        # ε_alg取得
        tEps_alg = self.mining_eps_alg(aTau)
        tEps_alg_size = len(tEps_alg.index)

        # q/ε_algでボンフェローニ補正を行う
        tPes = []
        for i in tEps_alg.index:
            tPes.append(tEps_alg_size*self.__pe(tEps_alg['Nep'][i],tEps_alg['Ne'][i]))
        tPv = pd.DataFrame(tPes,index=tEps_alg.index,columns=['p-value'])
        return pd.concat([tEps_alg,tPv],axis=1)

    def extract(self):
        """
        Mining and Multiple testing
        
        Returns
        -------
        tPtn_Pv : DataFrame  
            pattern and P-value
        """
        print("Step1,2を行う(Tarone的手順)")
        tTau = self.__get_min_tau()
        print("selcted tau:"+str(tTau))
        tPtn_Pv = self.__get_corrected_pvs(tTau)
        return tPtn_Pv

    def __kl(self,aP,aQ):
        """
        KLダイバージェンス計算
        """
        return aP*math.log(aP/aQ)+(1-aP)*math.log((1-aP)/(1-aQ))

    def __pe(self,aNep,aNe):
        """
        approximate P-value
        """
        tMue=aNep/aNe
        if tMue > self.mA:
            return math.exp(-aNe*self.__kl(tMue,self.mA))
        else:
            return 1

    #####仮 あとで別のファイルに分ける(とりあえずtestデータ)###
    # NeがaMin_sup以上のトランザクションを返す
    def mining_eps_alg(self, aMin_sup):
        return self.mDataset[self.mDataset['Ne']>= aMin_sup]





    