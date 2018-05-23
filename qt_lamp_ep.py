import math
import pandas as pd
import numpy as np

class Qt_Lamp_ep():
    def __init__(self,aA=0.7,aQ=0.05,aMax_min_sup=100,aMethod='BH',aDataset_carib=pd.DataFrame(),aDataset_main=pd.DataFrame()):
        """
        a in (0,1)

        Parameters
        -------
        aMax_min_sup : int
            min_supの最大値

        aDataset_main : DataFrame
            カラム名指定(['pattern','Nep','Ne'])
        aDataset_carib : DataFrame
            カラム名指定(['pattern','Nep','Ne'])
        """
        if(aA >= 1 or aA <= 0):
            print("Error:a=!(0~1)")
            exit()
        # if(aMethod =='BY' or aMethod =='BH'):
        #     print("Error:BH or BY")
        #     exit()
        self.mA = aA
        self.mMax_min_sup = aMax_min_sup
        self.mQ = aQ
        self.mDataset_carib = aDataset_carib
        self.mDataset_main = aDataset_main
        self.mMethod = aMethod

        
    def __get_min_tau(self):
        """
        τの一番小さいやつを得る
        """
  
        for i in range(1,self.mMax_min_sup):
            
            tEps_alg_carib = self.__mining_eps_alg(i,self.mDataset_carib)
            if(len(tEps_alg_carib.index)== 0):
                break
            # print("i="+str(i)+" "+str(self.mA**i)+"<="+str((self.mQ*self.__k(tEps_alg_carib))/(len(tEps_alg_carib.index)*self.__c(len(tEps_alg_carib.index))))+" k:"+str(self.__k(tEps_alg_carib))+" c_m:"+str(self.__c(len(tEps_alg_carib.index)))+" ceib_zize:"+str(len(tEps_alg_carib.index)))
            if self.mA**i <= (self.mQ*self.__k(tEps_alg_carib))/(len(tEps_alg_carib.index)*self.__c(len(tEps_alg_carib.index))):
                return i
        print("τ wasn't obtained")
        exit()

    def __get_corrected_pvs(self, aTau):
        """
        選ばれたτでパターン抽出を行い全てで多重検定をい,パターンとp値を返す
        """
        # ε_alg_main取得
        tEps_alg_main = self.__mining_eps_alg(aTau,self.mDataset_main)

        # q/ε_algでStep-up法を行う
        tK = self.__k(tEps_alg_main)
        print("selected k:"+str(tK))

        # p-valueを計算
        tPes = self.__pes(tEps_alg_main)
        tPv = pd.DataFrame(tPes,index=tEps_alg_main.index,columns=['p-value'])
        tEps_alg_main = pd.concat([tEps_alg_main,tPv],axis=1)

        # p-valueでソート
        tEps_alg_main = tEps_alg_main.sort_values(by='p-value') 
        
        if tK != 0:
            return tEps_alg_main[:tK]
        else:
            return pd.DataFrame(columns=tEps_alg_main.columns)



    def extract(self):
        """
        Mining and Multiple testing
        
        Returns
        -------
        tPtn_Pv : DataFrame  
            pattern and P-value
        """
        # print("Step1,2を行う(Tarone的手順)") 
        tTau = self.__get_min_tau()
        print("selcted tau:"+str(tTau))
        tPtn_Pv = self.__get_corrected_pvs(tTau)
        return tPtn_Pv

    def __kl(self,aP,aQ):
        """
        KLダイバージェンス計算
        """

        return aP*math.log(aP/aQ)+(1-aP)*math.log((1-aP)/(1-aQ))

    def __pes(self,aDataframe):
        """
        approximate P-value
        """

        tPes = []
        for i in aDataframe.index:
            tPes.append(self.__pe(aDataframe['Nep'][i],aDataframe['Ne'][i]))
        return tPes

    def __pe(self,aNep,aNe):
        tMue=aNep/aNe
        if tMue > self.mA:
            return math.exp(-aNe*self.__kl(tMue,self.mA))
        else:
            return 1

    def __c(self,aM):
        if(self.mMethod == 'BH'):
            return 1
        else:
            tSum = 0
            for i in range(1,aM+1):
                tSum = tSum + 1/i
            return tSum

    def __k(self,aDataframe):
        tPes = self.__pes(aDataframe)
        
        tPes.sort()
        tM=len(aDataframe.index)
        tC_m = self.__c(tM)  
        for i in range(1,tM+1):
            if(tPes[i-1] > (self.mQ*i)/(tC_m*tM)):                
                return i - 1        
        return tM

    # NeがaMin_sup以上のトランザクションを返す
    def __mining_eps_alg(self, aMin_sup,aDataframe):
        return aDataframe[aDataframe['Ne']>= aMin_sup]
