import tkinter as tk
from tkinter import ttk, filedialog
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf
from sklearn.preprocessing import PolynomialFeatures

def open_file():
    file = filedialog.askopenfile(mode='r', defaultextension='.csv',filetypes=[('Text Files', '*.csv'),('Text Files', '*.txt'),('Text Files', '*.dat')])
    if file:
        content = file.read()
        file.close()
        return content
      

      


def text_to_DataFrame(data):
    data=data.split("\n")
    data=[i.split() for i in data]
    try:
        df=pd.DataFrame(data).astype(float)
        df=df.dropna()
        return df
    except:
        df=pd.DataFrame(data[1:]).astype(float)
        df=df.dropna()
        return df

def Scaling_exponent_and_amplitude(data,scale="lin-lin"):
    '''This function calculates a scaling exponet and corresponding amplitude 
        using linear regression. It takes input data in the form of two columns. 
        For example for the gyration raduis it takes in log(N), log(Rg) or log(N), log(Rg^2)
        a simple scaling form of Rg^2=A*N^2n is assumed 
        
        data -- pd.DataFrame
        scale -- is data supplyed in lin-lin or log-log format. Lin-lin by definition
        output -- tuple of tuples where each inner tuple contains two numbers ((A,error_A),(nu,error_nu))'''
    data=str(data)
    try:
        if type(data)==str:
            data=text_to_DataFrame(data)
        if scale=="lin-lin":
            X=np.array(np.log(data.iloc[:,0])).reshape((-1,1))
            y=np.array(np.log(data.iloc[:,1]))
        elif scale=="log-log":
            X=np.array(data.iloc[:,0]).reshape((-1,1))
            y=np.array(data.iloc[:,1])
        else:
            message="Please check if submission is correct"
            return (np.nan,np.nan,(np.nan,np.nan),(np.nan,np.nan),message)
        X = sm.add_constant(X)
        score=0
        max_cutoff=int(0.9*y.size)
    
        for i in range(max_cutoff):
            results = sm.OLS(y[i:], X[i:]).fit()
            if results.rsquared>score:
                score=results.rsquared
                start=i
            else:
                break
        A,nu=results.params
        error_A,error_nu,=results.bse
        A=np.exp(A)
        error_A=A*error_A
        message=(f"Scaling exponent = {nu:.5f} ± {error_nu:.5f}\n Amplitude = {A:.5f} ± {error_A:.5f}\n"+
                 f"Calculation done starting from a data point {start}\n with accuracy of {results.rsquared:.5f}")
        return (start,results.rsquared,(A,error_A),(nu,error_nu),message)  
    except:
        message="Parameters cannot be calculated please check if submission is correct"
        return (np.nan,np.nan,(np.nan,np.nan),(np.nan,np.nan),message)

def Scaling_form_RG_chain(data,order='2'):
    '''In general Gyration radius can be presented in the scaling form 
        Rg^2=(A+B*N^(-Delta))*N^2nu
        Here Delta and nu are universal scaling exponents and A,B depend on both architecture and details of the model or chemistry
        This function is desighned to take in N and Rg^2 or R_g as input and return the best approximation for the parameters and their errorbars

        data -- pd.DataFrame
        order -- ether 1 or 2 is the power of thw observable by defalt it is 2
        returns -- ((A,error_A),(b,error_b))
    '''
    if type(data)==str:
        data=text_to_DataFrame(data)
    model=LinearRegression()
    X=np.array(np.log(data.iloc[:,0])).reshape((-1,1))
    if order=="1": 
        y=np.array(np.log(data.iloc[:,1]**2))
    else:
        y=np.array(np.log(data.iloc[:,1]))
    diff=0.1
    start=0
    max_cutoff=int(0.9*y.size)
    for i in range(max_cutoff):
        model.fit(X[i:],y[i:])
        if abs(model.coef_-1.175194)<diff:
            diff=abs(model.coef_-1.175194)
            start=i
        else:
            break
    coef=model.coef_
    try:
        while abs(coef-1.175194)>0.001:
            end=-1
            if coef-1.175194<0:
                model.fit(X[start:end],y[start:end])
                if abs(model.coef_-1.175194)<diff: 
                    diff=abs(model.coef_-1.175194)
                else:
                    end=end-1
                coef=model.coef_
            else:
                start=start+1
                model.fit(X[start:end],y[start:end])
                if abs(model.coef_-1.175194)<diff: 
                    diff=abs(model.coef_-1.175194)
                else:
                    start=start+1
                coef=model.coef_
            if diff-abs(coef-1.175194)<0.0001:
                print(diff)
                break
        X1=np.array(data.iloc[:,0][start:end]**(-0.528)).reshape((-1,1))
        y1=np.array(data.iloc[:,1][start:end]/data.iloc[:,0][start:end]**(1.175194))
        X1 = sm.add_constant(X1)
        results = sm.OLS(y1, X1).fit()
        A,b=results.params
        error_A,error_b=results.bse
        message=(f"Coeficients in scaling form:\n A = {A:.5f} ± {error_A:.5f}\n B = {b:.5f} ± {error_b:.5f}\n")
        return ((A,error_A),(b,error_b),message)  
    except:
        message="Parameters cannot be calculated please check if submission is correct"
        return ((np.nan,np.nan),(np.nan,np.nan),message)
def Rg_chain_prediction(N,sim_type=None,data=None):
    '''This function gives a value of Rg^2 of a chain given a number of monomers or molecular weight (N)
       NOTE eather type or data have to be provided and N is mandatory
       N -- a float number
       type -- a sting from a list (MC,MD) or None
       data -- pd.DataFrame
       returns Rg,error_Rg'''
    N=float(N)
    if type(data)==str:
        data=text_to_DataFrame(data)
    RG={"MC(build-in)":lambda n: 0.19514*(1-0.1125*n**(-0.528))*n**(1.175194),
        "MD(build-in)":lambda n:(0.26689-0.17305*n**(-0.528))*n**(1.175194)}
    error={"MC(build-in)":lambda n: RG["MC(build-in)"](n)*(2*np.log(n)*0.00007+(0.00004/0.19514)+abs(n**(-0.528)/((1-0.1125*n**(-0.528))))*0.0125
                                                           +abs((0.1125*n**(-0.528)*np.log(n))/((1-0.1125*n**(-0.528))))*0.012),
        "MD(build-in)":lambda n: RG["MD(build-in)"](n)*(2*np.log(n)*0.00007+(0.00180/(0.26689-0.17305*n**(-0.528)))+abs(n**(-0.528)/((0.26689-0.17305*n**(-0.528))))*0.03294
                                                           +abs((0.17305*n**(-0.528)*np.log(n))/((0.26689-0.17305*n**(-0.528))))*0.012)}
    if sim_type=="From my data": 
        A,b,_=Scaling_form_RG_chain(data)
        if A==np.nan and b==np.nan:
            message="Parameters cannot be calculated please check if submission is correct"
            return np.nan,np.nan,message
        A,error_A=A
        b,error_b=b
        Rg=(A+b*N**(-0.528))*N**(1.175194)
        error=Rg*(2*np.log(N)*0.00007+(error_A/(A+b*N**(-0.528)))+abs(N**(-0.528)/(A+b*N**(-0.528)))*error_b
                                                           +abs((0.17305*N**(-0.528)*np.log(N))/(A+b*N**(-0.528)))*0.012)
        message=f"Gyration radius squared = {Rg:.5f}± {error:.5f}.\n This result is based on user input data"
        return Rg,error,message
    if sim_type!="From my data":
        try:
            Rg=RG[sim_type](N)
            error=error[sim_type](N)
            
            if sim_type=="MC(build-in)":
                message=f"Gyration radius squared = {Rg:.5f}± {error:.5f}.\n This function was calculated using data from \n N.Clisby, PRL 104, 055702 (2010)"
            else:
                message=f"Gyration radius squared = {Rg:.5f}± {error:.5f}.\n This data was calculated from simulation results conducted for paper\n K. Haydukivska, V. Blavatska, and J. Paturej, PRE 108, 034502 (2023)"
            return Rg,error,message
        except:
            return np.nan,np.nan,"Input data provided incorrectly"
def Autocorrelation(data,mol_num):
    '''This function is desighned to provide a calculation of autocorreration function and fron it the relaxation time
    it uses three different aproaches for calculations
    input a time series to be correlated 
        
        data -- np.array

        tau -- itteration step float
        dump_step -- number of simulation steps between recording of the data
        mol_number -- number of moleules in the box
        time -- np.array of time values if avaliable


        return -- correlation time, max time of date'''
    mol_num=int(mol_num)
    if type(data)==str:
        data=text_to_DataFrame(data)
    try:
        time=np.array(data.iloc[:,0])
        rg_r=np.array(data.iloc[:,1:])
        corr=np.zeros(rg_r.shape)
        for i in range(mol_num):
            corr[:,i]=acf(rg_r[:,i],nlags=rg_r.shape[0])
        cor=np.mean(corr,axis=1)
        rel_time=np.trapz(cor[cor>0],time[cor>0])
        return rel_time,time[-1],f"Relaxation time = {rel_time}\n Full time = {time[-1]}"
    except:
        return np.nan,np.nan,"Please check if the data is correct"

def g1_g3(data,step=100,minsize=5000,accuracy=0.05):
    '''This function is designed to look for an exponent of 1 as a final regime of dynamic fuctions g1 and g3
    Input
    data -- a pd.DataFrame that contains two columns one for time second for g1 or g3 in linear scale
    step -- a size of the step in the loop by definition is 100
    minsize --- minimum number of points that are consideered or the fit
    accuracy -- a difference between the received exponent and 1'''
    step=int(step)
    minsize=int(minsize)
    accuracy=float(accuracy)
    if type(data)==str:
        data=text_to_DataFrame(data)
    try:
        model=LinearRegression()
        X=np.array(np.log(data.iloc[1:,0])).reshape((-1,1))
        y=np.array(np.log(data.iloc[1:,1]))
        diff=0.5
        start=0
        for i in range(1,y.size-minsize,step):
            model.fit(X[i:],y[i:])
            if abs(model.coef_-1.0)<diff: 
                if abs(model.coef_-1.0)>accuracy:
                    diff=abs(model.coef_-1.0)
                    start=i
                    score=model.score(X[i:],y[i:])
                else:
                    break
        score=model.score(X[start:],y[start:])
        return start,model.coef_,score,f"Starting from point {start} \n linear approximation gives an exponent {model.coef_}\n with model fitting accuracy of {score}"
    except:
        return np.nan,np.nan,"Please check if the data is correct"



def Is_there_scaling(data,step=1):
    """This function takes in data and returns approximated coefficients for whether:
            --- data.col_2=A*data.col_1^alpha 
            or
            --- data.col_2=A*exp(alpha*data.col_1)
            or 
            --- neither of those are true

        returns fitting parameters, range of column-1 for which the fitting was done or none if function was not found
    """
    step=int(step)
    if type(data)==str:
        data=text_to_DataFrame(data)
    X=np.array(np.log(data.iloc[:,0])).reshape((-1,1))
    y=np.array(np.log(data.iloc[:,1]))
    X = sm.add_constant(X)
    max_cutoff=int(0.9*y.size)
    power_law=False
    try:
        score=0
        for i in range(0,max_cutoff,step):
            results = sm.OLS(y[i:], X[i:]).fit()
            if results.rsquared>score:
                score=results.rsquared
                start=i
            else:
                break
        if len(X[start:])>4:
            power_law=True
            A,nu=results.params
            error_A,error_nu=results.bse
            A=np.exp(A)
            error_A=A*error_A
            message=f"Calculations conducted successfully. A power law found.\n Amplitude = {A:.5f}± {error_A:.5f} \n Exponent = {nu:.5f}± {error_nu:.5f} \n Starting from point {start}\n calculated with accuracy {score}"
        else:  
            A=np.nan;error_A=np.nan;nu=np.nan;error_nu=np.nan
            message="There is not enough data to conclusivly state power law dependence"
        if power_law==True: 
            return ((A,error_A),(nu,error_nu),message)
        X=np.array(data.iloc[:,0]).reshape((-1,1))
        X = sm.add_constant(X)
        score=0
        for i in range(0,max_cutoff,step):
            results = sm.OLS(y[i:], X[i:]).fit()
            if results.rsquared>score:
                score=results.rsquared
                start=i
            else:
                break
        if len(X[start:])>4:
            A,nu=results.params
            error_A,error_nu,=results.bse
            A=np.exp(A)
            error_A=A*error_A
            message=f"Calculations conducted successfully. An exponential dependence found.\n Amplitude = {A:.5f}± {error_A:.5f} \n Exponent parameter = {nu:.5f}± {error_nu:.5f} \n Starting from point {start}\n calculated with accuracy {score}"
        else:  
            A=np.nan;error_A=np.nan;nu=np.nan;error_nu=np.nan
            message="There is not enough data to conclusivly state exponential dependence"
        return ((A,error_A),(nu,error_nu),message)
    except:
        message="Parameters cannot be calculated please check if submission is correct"
        return ((np.nan,np.nan),(np.nan,np.nan),message)

def finit_size_scaling(data):
    '''The function is designed to calculate a finit size scaling of size ratios and 
        aspherisity

        data -- pd.DataFrame
        returns -- ((A,error_A),(b,error_b))
    '''
    if type(data)==str:
        data=text_to_DataFrame(data)
    model=LinearRegression()
    poly = PolynomialFeatures(degree=2, include_bias=False)
    X=np.array(data.iloc[:,0]**(-0.528))
    poly_features = poly.fit_transform(X.reshape(-1, 1))
    y=np.array(data.iloc[:,1])
    score=0.0
    start=0
    max_cutoff=int(0.9*y.size)
    for i in range(max_cutoff):
        model.fit(poly_features[i:],y[i:])
        if model.score(poly_features[i:],y[i:])>score:
            score=model.score(poly_features[i:],y[i:])
            start=i
        else:
            continue
    model.fit(poly_features[start:],y[start:])
    X1=np.hstack((np.array(data.iloc[:,0][start:]**(-0.528)).reshape((-1,1)),np.array(data.iloc[:,0][start:]**(-0.528*2)).reshape((-1,1))))
    y1=np.array(data.iloc[:,1][start:])
    try:
        X1 = sm.add_constant(X1)
        results = sm.OLS(y1, X1).fit()
        A,b,c=results.params
        error_A,error_b,error_c=results.bse
        return (start,(A,error_A),(b,error_b),f"For a infinitly long chain the parameter is = {A:.5f}± {error_A:.5f}" )
    except:
        message="Parameters cannot be calculated please check if submission is correct"
        return (np.nan,(np.nan,np.nan),(np.nan,np.nan),message)