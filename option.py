# -*- coding: cp936 -*-
#python version: 3.7.4
#本程序主要用于计算期权定价（基于black-scholes方程）
#计算欧式看涨/跌期权和美式看涨/跌期权四种期权的定价
#对于欧式期权采取真实解和二叉树方法求近似解
#对于美式期权采取二叉树方法和最小二乘蒙特卡洛方法进行求解

from math import log,sqrt,exp
from scipy.stats import norm
import numpy as np
from numpy.linalg import inv
from tkinter import *
from tkinter.messagebox import *
from tkinter import ttk

class option:
    def __init__(self,cur_price, exe_price,time,vol = 0.2,arbi=0.1,n = 6,m = 500):
        self.current_price = cur_price #股票现价
        self.executive_price = exe_price #行权价格
        self.time = time  #期权时效
        self.volatility = vol  #股价波动率
        self.riskless_arbitrage = arbi   #无风险利率
        self.option = {}  #算得相应期权的定价
        self.n = n  #将总体时间time分成n小段进行二叉树算法的计算
        self.m = m  #在最小二乘蒙特卡洛中表示结果是经历了m次随机取平均的结果
    
    def Information(self):
        print('Current Price:',self.current_price)
        print('Executive Price:',self.executive_price)
        print('Time:',self.time)
        print('Volatility:',self.volatility)
        print('Riskless Arbitrage:',self.riskless_arbitrage)
        print('n:',self.n)
        print('m:',self.m)
        
    #对欧式期权计算精确解
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

    #对欧式期权用二叉树方法求解
    #这里考虑T以年为单位（为了和利率单位一致）,分成n个时间段进行计算
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

    #对美式期权用二叉树方法求解
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

    #一个工具函数，用于最小二乘蒙塔卡罗中求解最小二乘的解
    def pre(self,S,m,i):
        A = np.ones([m,3])
        A[:,0] = 1
        A[:,1] = S[:,i]
        A[:,2] = S[:,i]*S[:,i]
        y = S[:,i+1]
        z = np.dot(np.dot(inv(np.dot(A.T,A)),A.T),y)
        return np.dot(A,z)

    #在美式期权中用最小二乘蒙特卡洛求解
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
    #这样方便下面进行调用
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

#程序的类
class MyApp:
    dictt = {'欧式看涨期权 European Call Option':'EC',
                 '欧式看跌期权 European Put Option':'EP',
                 '美式看涨期权 American Call Option':'AC',
                 '美式看跌期权 American Put Option':'AP',
                 'Black-Scholes精确解':'P',
                 '二叉树模拟':'BT',
                 '最小二乘蒙特卡洛模拟':'LSM'}
    def __init__(self):
        self.root = Tk()
        self.root.title('期权定价计算器')
        self.root.geometry('680x170+200+200')
        
        #初始化参数
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
        
        #基本输入框
        self.cpl = Label(self.root,text = '股票现价:',width = 16,anchor=W,justify="left")
        self.cpl.grid(sticky = W)
        self.epl = Label(self.root,text = '行权价格:')
        self.epl.grid(sticky = W)
        self.rl = Label(self.root,text = '无风险利率:')
        self.rl.grid(sticky = W)
        self.cp = Entry(self.root,textvariable = self.cur_price,width = 30)
        self.cp.grid(row = 0, column = 1,sticky = W)
        self.ep = Entry(self.root,textvariable = self.exe_price,width = 30)
        self.ep.grid(row = 1, column = 1,sticky = W)
        self.r = Entry(self.root,textvariable = self.riskless_arbitrage,width = 30)
        self.r.grid(row = 2, column = 1,sticky = W)
        self.vl = Label(self.root,text = '股票波动率:',width = 16,anchor=W,justify="left")
        self.vl.grid(sticky = W,row = 0, column = 2)
        self.tl = Label(self.root,text = '期权时效（年）:')
        self.tl.grid(sticky = W,row = 1, column = 2)
        self.typl = Label(self.root,text = '期权类型:')
        self.typl.grid(sticky = W,row = 2, column = 2)
        self.v = Entry(self.root,textvariable = self.votatility,width = 30)
        self.v.grid(row = 0, column = 3,sticky = W)
        self.t = Entry(self.root,textvariable = self.time,width = 30)
        self.t.grid(row = 1, column = 3,sticky = W)
        #计算期权类型下拉框
        self.typ = ttk.Combobox(self.root, width=27,state = 'readonly')
        self.typ.grid(row = 2, column = 3, sticky = W)
        self.typ['value'] = ('欧式看涨期权 European Call Option',
                              '欧式看跌期权 European Put Option',
                              '美式看涨期权 American Call Option',
                              '美式看跌期权 American Put Option')
        self.typ.current(0)
        #计算方法下拉框（和期权类型联动，因为不同期权的计算方式不同）
        self.methodl = Label(self.root,text = '计算方式:')
        self.methodl.grid(sticky = W,row = 3, column = 0)
        self.method = ttk.Combobox(self.root,width = 27, state = 'readonly')
        #先初始化计算方法的框图
        if self.typ.get() == '欧式看涨期权 European Call Option' or self.typ.get() == '欧式看跌期权 European Put Option':
            self.method['value'] = ('Black-Scholes精确解','二叉树模拟')
        else:
            self.method['value'] = ('二叉树模拟','最小二乘蒙特卡洛模拟')
        self.method.current(0)
        #随着type变化而自动变化
        self.typ.bind("<<ComboboxSelected>>",self.changemethod)
        self.method.grid(row = 3, column = 1, sticky = W)
        
        #按照方法选择是否出现的框图
        self.nl = Label(self.root,text = '二叉树分割时间数:')
        self.n = Entry(self.root,textvariable = self.nv,width = 30)
        if self.method.get() != 'Black-Scholes精确解':
            self.nl.grid(row = 4,column = 0,sticky = W)
            self.n.grid(row = 4,column = 1,sticky = W)
        self.ml = Label(self.root,text = '蒙特卡洛重复次数:')
        self.m = Entry(self.root,textvariable = self.mv,width = 30)
        if self.method.get() == '最小二乘蒙特卡洛模拟':
            self.ml.grid(row = 4,column = 2,sticky = W)
            self.m.grid(row = 4,column = 3,sticky = W)
        self.method.bind("<<ComboboxSelected>>",self.changemn)

        #计算按键
        self.cal = Button(self.root,text = 'Calculate',command = self.Calculate)
        self.cal.grid(row = 5,column = 3)
        self.res = Label(self.root,text = '结果:')
        self.res.grid(row = 5,column = 0,sticky =W)
        self.resv = Entry(self.root,textvariable = self.result,state = 'readonly',width  =30)
        self.resv.grid(row = 5, column = 1, sticky = W)

        #退出提醒
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
        
    #计算方法的下拉框和对应参数随着选择期权类型改变
    def changemethod(self,*arg):
        if self.typ.get() == '欧式看涨期权 European Call Option' or self.typ.get() == '欧式看跌期权 European Put Option':
            self.method['value'] = ('Black-Scholes精确解','二叉树模拟')
        else:
            self.method['value'] = ('二叉树模拟','最小二乘蒙特卡洛模拟')
        self.method.current(0)
        
        if self.method.get() != 'Black-Scholes精确解':
            self.nl.grid(row = 4,column = 0,sticky = W)
            self.n.grid(row = 4,column = 1,sticky = W)
        else:
            self.nl.grid_forget()
            self.n.grid_forget()
        if self.method.get() == '最小二乘蒙特卡洛模拟':
            self.ml.grid(row = 4,column = 2,sticky = W)
            self.m.grid(row = 4,column = 3,sticky = W)
        else:
            self.ml.grid_forget()
            self.m.grid_forget()
    
    #对应计算方法的参数框格随着计算方法而改变
    def changemn(self,*arg):
        if self.method.get() != 'Black-Scholes精确解':
            self.nl.grid(row = 4,column = 0,sticky = W)
            self.n.grid(row = 4,column = 1,sticky = W)
        else:
            self.nl.grid_forget()
            self.n.grid_forget()
        if self.method.get() == '最小二乘蒙特卡洛模拟':
            self.ml.grid(row = 4,column = 2,sticky = W)
            self.m.grid(row = 4,column = 3,sticky = W)
        else:
            self.ml.grid_forget()
            self.m.grid_forget()
            
    #调用option类计算期权
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
