# -*- coding: cp936 -*-
#python version: 3.7.4
#��������Ҫ���ڼ�����Ȩ���ۣ�����black-scholes���̣�
#����ŷʽ����/����Ȩ����ʽ����/����Ȩ������Ȩ�Ķ���
#����ŷʽ��Ȩ��ȡ��ʵ��Ͷ�������������ƽ�
#������ʽ��Ȩ��ȡ��������������С�������ؿ��巽���������

from math import log,sqrt,exp
from scipy.stats import norm
import numpy as np
from numpy.linalg import inv
from tkinter import *
from tkinter.messagebox import *
from tkinter import ttk

class option:
    def __init__(self,cur_price, exe_price,time,vol = 0.2,arbi=0.1,n = 6,m = 500):
        self.current_price = cur_price #��Ʊ�ּ�
        self.executive_price = exe_price #��Ȩ�۸�
        self.time = time  #��ȨʱЧ
        self.volatility = vol  #�ɼ۲�����
        self.riskless_arbitrage = arbi   #�޷�������
        self.option = {}  #�����Ӧ��Ȩ�Ķ���
        self.n = n  #������ʱ��time�ֳ�nС�ν��ж������㷨�ļ���
        self.m = m  #����С�������ؿ����б�ʾ����Ǿ�����m�����ȡƽ���Ľ��
    
    def Information(self):
        print('Current Price:',self.current_price)
        print('Executive Price:',self.executive_price)
        print('Time:',self.time)
        print('Volatility:',self.volatility)
        print('Riskless Arbitrage:',self.riskless_arbitrage)
        print('n:',self.n)
        print('m:',self.m)
        
    #��ŷʽ��Ȩ���㾫ȷ��
    def EP(self,typ = 'C'):
        S = self.current_price
        X = self.executive_price
        r = self.riskless_arbitrage
        v = self.volatility
        T = self.time
        d1 = (log(S/X)+(r+0.5*v*v)*T)/v/sqrt(T)
        d2 = d1 - v*sqrt(T)
        if typ == 'C':
            return S*norm.cdf(d1)-exp(-r*T)*X*norm.cdf(d2)
        if typ == 'P':
            return -S*norm.cdf(-d1)+exp(-r*T)*X*norm.cdf(-d2)

    #��ŷʽ��Ȩ�ö������������
    #���￼��T����Ϊ��λ��Ϊ�˺����ʵ�λһ�£�,�ֳ�n��ʱ��ν��м���
    def EBT(self,typ = 'C'):
        v = self.volatility
        r = self.riskless_arbitrage
        t = self.time
        S0 = self.current_price
        X = self.executive_price
        n = self.n
        dt = t/n
        u = exp(v*sqrt(dt))
        d = 1/u
        p = (exp(r*dt)-d)/(u-d)
        S = list()
        V = list()
        for i in range(n+1):
            S.append(S0*u**i*d**(n-i))
            if typ == 'P':
                V.append(max(-S[i]+X,0))
            if typ == 'C':
                V.append(max(S[i]-X,0))
        for i in range(n):
            Vn = []
            for j in range(len(V)-1):
                Vn.append(exp(-r*dt)*((1-p)*V[j] + p*V[j+1]))
            V = Vn
        return V[0]

    #����ʽ��Ȩ�ö������������
    def ABT(self,typ = 'C'):
        v = self.volatility
        r = self.riskless_arbitrage
        t = self.time
        S0 = self.current_price
        X = self.executive_price
        n = self.n
        dt = t/n
        u = exp(v*sqrt(dt))
        d = 1/u
        p = (exp(r*dt)-d)/(u-d)
        S = list()
        V = list()
        for i in range(n+1):
            S.append([])
            for j in range(i+1):
                S[i].append(S0*u**j*d**(i-j))
        for i in range(n+1):
            if typ == 'C':
                V.append(max(S[n][i] - X,0))
            if typ == 'P':
                V.append(max(- S[n][i] + X,0))
        for i in range(n):
            Vn = list()
            for j in range(len(V)-1):
                if typ == 'C':
                    St = max(S[n-1-i][j]-X,0)
                if typ == 'P':
                    St = max(-S[n-1-i][j]+X,0)
                no = exp(-r*dt)*((1-p)*V[j]+p*V[j+1])
                Vn.append(max(St,no))
            V = Vn
        return V[0]

    #һ�����ߺ�����������С�������������������С���˵Ľ�
    def pre(self,S,m,i):
        A = np.ones([m,3])
        A[:,0] = 1
        A[:,1] = S[:,i]
        A[:,2] = S[:,i]*S[:,i]
        y = S[:,i+1]
        z = np.dot(np.dot(inv(np.dot(A.T,A)),A.T),y)
        return np.dot(A,z)

    #����ʽ��Ȩ������С�������ؿ������
    def ALSM(self,typ = 'C'):
        v = self.volatility
        r = self.riskless_arbitrage
        t = self.time
        S0 = self.current_price
        X = self.executive_price
        n = self.n
        m = self.m
        dt = t/n
        u = exp(v*sqrt(dt))
        d = 1/u
        p = (exp(r*dt)-d)/(u-d)
        
        S = np.zeros([m,n+1])
        S[:,0] = S0
        for i in range(m):
            for j in range(1,n+1):
                rand = np.random.randn(1)[0]
                S[i][j] = S[i][j-1]*exp((r - v**2*0.5)*dt+v*sqrt(dt)*rand)
        A = np.zeros([m,n+1])
        A[:,0] = S0
        A[:,1] = S[:,1]
        for i in range(2,n+1):
            A[:,i] = self.pre(S,m,i-1)
        if typ == 'C':
            A = A-X
            S = S-X
        if typ == 'P':
            A = X-A
            S = X-S
        A = np.array([[max(A[i][j],0) for j in range(A.shape[1])] for i in range(A.shape[0])])
        S = np.array([[max(S[i][j],0) for j in range(S.shape[1])] for i in range(S.shape[0])])
        
        summ = 0
        for i in range(m):
            k = 0
            for j in range(n):
                if A[i][j+1]<=S[i][j] and k==0:
                    k = j
            if k ==0:
                k = i+1
            summ += exp(-r*dt*(k-1))*S[i][j]
        summ = summ/m
        return summ
    
    ###=====================================================
    #��������������е���
    def ECP(self):
        return self.EP('C')
    def EPP(self):
        return self.EP('P')
    def ECBT(self):
        return self.EBT('C')
    def EPBT(self):
        return self.EBT('P')
    def ACBT(self):
        return self.EBT('C')
    def APBT(self):
        return self.EBT('P')
    def ACLSM(self):
        return self.ALSM('C')
    def APLSM(self):
        return self.ALSM('P')
    ###=======================================================
    def Calculate_option(self,typ = 'EC',method = 'P'):
        switch = {'ECP':self.ECP,
                  'EPP':self.EPP,
                  'ECBT':self.ECBT,
                  'EPBT':self.EPBT,
                  'ACBT':self.ACBT,
                  'APBT':self.APBT,
                  'ACLSM':self.ACLSM,
                  'APLSM':self.APLSM}
        return switch[typ+method]()

#�������
class MyApp:
    dictt = {'ŷʽ������Ȩ European Call Option':'EC',
                 'ŷʽ������Ȩ European Put Option':'EP',
                 '��ʽ������Ȩ American Call Option':'AC',
                 '��ʽ������Ȩ American Put Option':'AP',
                 'Black-Scholes��ȷ��':'P',
                 '������ģ��':'BT',
                 '��С�������ؿ���ģ��':'LSM'}
    def __init__(self):
        self.root = Tk()
        self.root.title('��Ȩ���ۼ�����')
        self.root.geometry('680x170+200+200')
        
        #��ʼ������
        self.cur_price = StringVar()
        self.cur_price.set('50.00')
        self.exe_price = StringVar()
        self.exe_price.set('50.00')
        self.time = StringVar()
        self.time.set('0.5000')
        self.riskless_arbitrage = StringVar()
        self.riskless_arbitrage.set('0.10')
        self.votatility = StringVar()
        self.votatility.set('0.20')
        self.nv = StringVar()
        self.nv.set('6')
        self.mv = StringVar()
        self.mv.set('1000')
        self.result = StringVar()
        self.result.set('')
        
        #���������
        self.cpl = Label(self.root,text = '��Ʊ�ּ�:',width = 16,anchor=W,justify="left")
        self.cpl.grid(sticky = W)
        self.epl = Label(self.root,text = '��Ȩ�۸�:')
        self.epl.grid(sticky = W)
        self.rl = Label(self.root,text = '�޷�������:')
        self.rl.grid(sticky = W)
        self.cp = Entry(self.root,textvariable = self.cur_price,width = 30)
        self.cp.grid(row = 0, column = 1,sticky = W)
        self.ep = Entry(self.root,textvariable = self.exe_price,width = 30)
        self.ep.grid(row = 1, column = 1,sticky = W)
        self.r = Entry(self.root,textvariable = self.riskless_arbitrage,width = 30)
        self.r.grid(row = 2, column = 1,sticky = W)
        self.vl = Label(self.root,text = '��Ʊ������:',width = 16,anchor=W,justify="left")
        self.vl.grid(sticky = W,row = 0, column = 2)
        self.tl = Label(self.root,text = '��ȨʱЧ���꣩:')
        self.tl.grid(sticky = W,row = 1, column = 2)
        self.typl = Label(self.root,text = '��Ȩ����:')
        self.typl.grid(sticky = W,row = 2, column = 2)
        self.v = Entry(self.root,textvariable = self.votatility,width = 30)
        self.v.grid(row = 0, column = 3,sticky = W)
        self.t = Entry(self.root,textvariable = self.time,width = 30)
        self.t.grid(row = 1, column = 3,sticky = W)
        #������Ȩ����������
        self.typ = ttk.Combobox(self.root, width=27,state = 'readonly')
        self.typ.grid(row = 2, column = 3, sticky = W)
        self.typ['value'] = ('ŷʽ������Ȩ European Call Option',
                              'ŷʽ������Ȩ European Put Option',
                              '��ʽ������Ȩ American Call Option',
                              '��ʽ������Ȩ American Put Option')
        self.typ.current(0)
        #���㷽�������򣨺���Ȩ������������Ϊ��ͬ��Ȩ�ļ��㷽ʽ��ͬ��
        self.methodl = Label(self.root,text = '���㷽ʽ:')
        self.methodl.grid(sticky = W,row = 3, column = 0)
        self.method = ttk.Combobox(self.root,width = 27, state = 'readonly')
        #�ȳ�ʼ�����㷽���Ŀ�ͼ
        if self.typ.get() == 'ŷʽ������Ȩ European Call Option' or self.typ.get() == 'ŷʽ������Ȩ European Put Option':
            self.method['value'] = ('Black-Scholes��ȷ��','������ģ��')
        else:
            self.method['value'] = ('������ģ��','��С�������ؿ���ģ��')
        self.method.current(0)
        #����type�仯���Զ��仯
        self.typ.bind("<<ComboboxSelected>>",self.changemethod)
        self.method.grid(row = 3, column = 1, sticky = W)
        
        #���շ���ѡ���Ƿ���ֵĿ�ͼ
        self.nl = Label(self.root,text = '�������ָ�ʱ����:')
        self.n = Entry(self.root,textvariable = self.nv,width = 30)
        if self.method.get() != 'Black-Scholes��ȷ��':
            self.nl.grid(row = 4,column = 0,sticky = W)
            self.n.grid(row = 4,column = 1,sticky = W)
        self.ml = Label(self.root,text = '���ؿ����ظ�����:')
        self.m = Entry(self.root,textvariable = self.mv,width = 30)
        if self.method.get() == '��С�������ؿ���ģ��':
            self.ml.grid(row = 4,column = 2,sticky = W)
            self.m.grid(row = 4,column = 3,sticky = W)
        self.method.bind("<<ComboboxSelected>>",self.changemn)

        #���㰴��
        self.cal = Button(self.root,text = 'Calculate',command = self.Calculate)
        self.cal.grid(row = 5,column = 3)
        self.res = Label(self.root,text = '���:')
        self.res.grid(row = 5,column = 0,sticky =W)
        self.resv = Entry(self.root,textvariable = self.result,state = 'readonly',width  =30)
        self.resv.grid(row = 5, column = 1, sticky = W)

        #�˳�����
        self.root.protocol('WM_DELETE_WINDOW',self.callback)
        self.root.update()
        self.root.mainloop()
        
    def callback(self):
        if askokcancel('Quit','Do you really want to quit?'):
            self.root.destroy()
    def printmethod(self,*arg):
        print(self.method.get())
    def printtype(self,*arg):
        print(self.typ.get())
        
    #���㷽����������Ͷ�Ӧ��������ѡ����Ȩ���͸ı�
    def changemethod(self,*arg):
        if self.typ.get() == 'ŷʽ������Ȩ European Call Option' or self.typ.get() == 'ŷʽ������Ȩ European Put Option':
            self.method['value'] = ('Black-Scholes��ȷ��','������ģ��')
        else:
            self.method['value'] = ('������ģ��','��С�������ؿ���ģ��')
        self.method.current(0)
        
        if self.method.get() != 'Black-Scholes��ȷ��':
            self.nl.grid(row = 4,column = 0,sticky = W)
            self.n.grid(row = 4,column = 1,sticky = W)
        else:
            self.nl.grid_forget()
            self.n.grid_forget()
        if self.method.get() == '��С�������ؿ���ģ��':
            self.ml.grid(row = 4,column = 2,sticky = W)
            self.m.grid(row = 4,column = 3,sticky = W)
        else:
            self.ml.grid_forget()
            self.m.grid_forget()
    
    #��Ӧ���㷽���Ĳ���������ż��㷽�����ı�
    def changemn(self,*arg):
        if self.method.get() != 'Black-Scholes��ȷ��':
            self.nl.grid(row = 4,column = 0,sticky = W)
            self.n.grid(row = 4,column = 1,sticky = W)
        else:
            self.nl.grid_forget()
            self.n.grid_forget()
        if self.method.get() == '��С�������ؿ���ģ��':
            self.ml.grid(row = 4,column = 2,sticky = W)
            self.m.grid(row = 4,column = 3,sticky = W)
        else:
            self.ml.grid_forget()
            self.m.grid_forget()
            
    #����option�������Ȩ
    def Calculate(self):
        self.option = option(eval(self.cur_price.get()),
                             eval(self.exe_price.get()),
                             eval(self.time.get()),
                             eval(self.votatility.get()),
                             eval(self.riskless_arbitrage.get()),
                             eval(self.nv.get()),
                             eval(self.mv.get()))
        self.result.set(str(round(self.option.Calculate_option(self.dictt[self.typ.get()],self.dictt[self.method.get()]),2)))
        
if __name__ == '__main__':
    MyApp()
