import scipy.signal as ssig
import math
import numpy as np

class QRS_detector:
    
    def __init__(self,name):
        self.name=name
    
    def load_data(self): # load data
        file=open(self.name,'r')
        lines=file.readlines()
        time,ecg=[],[]
        for i in lines:
            time.append(float(i.split(" ")[2]))
            ecg.append(float(i.split(" ")[-1]))
        return time, ecg
    
    def high_Freq_response(self,data,fs): # Highpass_filter + freq_response
        b,a= ssig.butter(6,10,btype='high',fs=fs)
        w,h=ssig.freqz(b,a,fs)
        w_hz=w*fs/(2*math.pi)
        high_data=ssig.filtfilt(b,a,data)
        return w_hz, abs(h), high_data
    
    def low_Freq_response(self,data,fs): # Lowpass_filter + freq_response
        b,a=ssig.butter(5,50,btype='low', fs=fs)
        w,h=ssig.freqz(b,a,fs)
        w_hz=w*fs/(2*math.pi)
        low_data=ssig.filtfilt(b,a,data)
        return w_hz, abs(h), low_data
    
    def Differentiation(self,data,h): # Differentiation 
        data=np.array(data)
        return np.r_[0, np.diff(data)]/h
    
    def Moving_Average(self,data): # moving average filter 
        Moving_result=[0]*8
        for i in range(8,len(data)):
            Moving_result.append(sum(data[i-8:i])/8)
        return Moving_result
    
    def peak_detector_class(self,data,fs): # peak detector 
        '''
        예상 피크점 검출시, 뒤 세개의 점의 평균을 계산 후, 크기 여부 확인하여 피크 확정 검출 
        '''
        threshold=0.6*max(data)
        pre_peak,ori_peak=[0]*len(data),[0]*len(data)
        for i in range(len(data)):
            if data[i]>data[i-1] and data[i]> data[i+1] and threshold< data[i]:
                if data[i]> sum([data[i+1],data[i+2],data[i+3]])/3:
                   pre_peak[i],ori_peak[i]=1,data[i]
        return ori_peak,pre_peak
    
    def bpm_calculator(self,data,time,fs): # bpm 
        number,dsd,bpm=[],[],[]
        for i in data:
            if i !=0:
                number.append(time[data.index(i)])
        for tm in range(0,len(time),400):
            ex_time=time[tm:tm+400]
            for d in number:
                if d in ex_time:
                    dsd.append(d)
            rri_interval=[dsd[i+1]-dsd[i] for i in range(len(dsd)-1)]
            dsd=[]
            bpm.append( 60 /((sum(rri_interval))/len(rri_interval)))
        bpm_time=np.linspace(0,8,len(bpm)) 
        return bpm_time, bpm
        
    def show_result(self): 
        time,ecg=QRS_detector.load_data(self)
        fs=int(1/(time[1]-time[0]))
        w_hz,h_fr,high_ecg=QRS_detector.high_Freq_response(self,ecg,fs)
        l_hz,l_fr,low_ecg=QRS_detector.low_Freq_response(self,high_ecg,fs) # result of Bandpassfilter
        diff=QRS_detector.Differentiation(self,low_ecg,1/fs) # result of Differentiation 
        squaring_data=[i**2 for i in low_ecg] # result of squaring 
        moving_result=QRS_detector.Moving_Average(self,squaring_data) # result of moving average 
        peak_detector,ori_peak=QRS_detector.peak_detector_class(self,moving_result,fs)
        bpm_time,bpm=QRS_detector.bpm_calculator(self,peak_detector,time,fs)
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,10))
        plt.subplot(4,1,1);plt.plot(time,diff,'black');plt.title("Result of Differentiation");plt.xlabel('Time(s)'); plt.ylabel('Amplitude(mV)');
        plt.subplot(4,1,2);plt.plot(time,squaring_data,'black');plt.title("Result of Squaring");plt.xlabel('Time(s)'); plt.ylabel('Amplitude(mV)');
        plt.subplot(4,1,3);plt.plot(time,moving_result,'black');plt.title("Result of Moving average");plt.xlabel('Time(s)'); plt.ylabel('Amplitude(mV)');
        plt.subplot(4,2,7);plt.plot(time,ecg,'black');plt.plot(time,ori_peak,'red');plt.title("Peak Detection");plt.xlabel('Time(s)'); plt.ylabel('Amplitude(mV)');
        plt.subplot(4,2,8);plt.bar(bpm_time,bpm, color="black");plt.title("Heart Rate(t=2s interval)");plt.xlabel('Time(s)'); plt.ylabel('Heart Rate(bpm)');
        plt.tight_layout();plt.show()

if __name__=="__main__":
    QRS_detector('ECG_data.txt').show_result()

    