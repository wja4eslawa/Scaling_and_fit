import tkinter as tk
from tkinter import ttk, filedialog
from Functions import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

def calculate(function,input,*output):
    """The button function depending on the input calculates a result and outputs it into separate window
    Takes in 
    function -- str which is a name of a function in a file Functions
    input -- a tuple of strings taken from input fields
    output -- outputs of the calculations they are st to corresponting variables in the interfase"""
    if function=="Scaling_exponent_and_amplitude":
        data,scale=input
        returns=Scaling_exponent_and_amplitude(data,str(scale))
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
        
    elif function=="Scaling_form_RG_chain":
        data,order=input
        returns=Scaling_form_RG_chain(data,order)
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
    elif function=="Rg_chain_prediction":
        N,type,data=input
        returns=Rg_chain_prediction(N,type,data)
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
    elif function=="finit_size_scaling":
        data=input
        returns=finit_size_scaling(data)
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
    elif function=="Is_there_scaling":
        data,step=input
        returns=Is_there_scaling(data,step)
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
    elif function=="g1_g3":
        data,step,minsize,accuracy=input
        returns=g1_g3(data,step,minsize,accuracy)
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
    else: 
        data,mol_num=input
        returns=Autocorrelation(data,mol_num)
        newWindow = tk.Toplevel(wind)
        newWindow.title("Results of the calculation")
        newWindow.geometry("400x100")
        tk.Label(newWindow,text=returns[-1]).pack()
    for i in range(len(returns)):
        output[i].set(returns[i])


def display_picture(function,start):
      """This function is plotting the results of the calculation"""
      newWindow = tk.Toplevel(wind)
      newWindow.title("Fit to the data")
      newWindow.geometry("600x600")
      data=tk.StringVar(value="data")

      file = filedialog.askopenfile(mode='r', defaultextension='.csv',filetypes=[('Text Files', '*.csv'),('Text Files', '*.txt'),('Text Files', '*.dat')])
      if file:
            content = file.read()
            file.close()
            data.set(content)
      fig,ax=plt.subplots()
      canvas=FigureCanvasTkAgg(fig,master=newWindow)
      ##input data to the figure
      input_data=text_to_DataFrame(data.get())
      start=int(start.get())
      if function=="Scaling_exponent_and_amplitude":
            x=np.log(np.array(input_data[0]))
            y=np.log(np.array(input_data[1]))
            model=LinearRegression()
            model.fit(x[start:].reshape((-1,1)),y[start:])
            if y.size<1000:
                ax.scatter(x,y,color="black",label="data",s=5)
            else:
               ax.plot(x,y,color="black",label="data") 
            ax.plot(x[start:],model.predict(x[start:].reshape((-1,1))),color='red',label="fit")
            plt.legend()
            plt.xlabel("Natural log of number of monomers")
            plt.ylabel("Natural log of size characteristic")
      elif function=="finit_size_scaling":
            x=np.array(input_data[0])**(-0.528)
            y=np.array(input_data[1])
            model=LinearRegression()
            poly = PolynomialFeatures(degree=2, include_bias=False)
            poly_features = poly.fit_transform(x.reshape(-1, 1))
            model.fit(poly_features[start:],y[start:])
            ax.scatter(x,y,color="black",label="data",s=5)
            ax.plot(x[start:],model.predict(poly_features[start:]),color='red',label="fit")
            plt.legend()
            plt.xlabel("$N^{-\\Delta}$")
            plt.ylabel("Universal characteristic")
      else:
            x=np.log(np.array(input_data[0])[1:])
            y=np.log(np.array(input_data[1])[1:])
            model=LinearRegression()
            model.fit(x[start:].reshape((-1,1)),y[start:])
            if y.size<1000:
                  ax.scatter(x,y,color="black",label="data",s=5)
            else:
                  ax.plot(x,y,color="black",label="data") 
            ax.plot(x[start:],model.predict(x[start:].reshape((-1,1))),color='red',label="fit")
            plt.legend()
            plt.xlabel("Natural log of time")
            plt.ylabel("Natural log of dynamic characteristic")
      canvas.draw()
      toolbar=NavigationToolbar2Tk(canvas,newWindow,pack_toolbar=False)
      toolbar.update()
      toolbar.pack(anchor="w",fill=tk.X)

      canvas.get_tk_widget().pack()

def data_preparation():
      """The function provides instructions on the file preparations that apperas in the separate window"""
      newWindow = tk.Toplevel(wind)
      newWindow.title("Results of the calculation")
      newWindow.geometry("400x400")
      message="""The data for the functions has to be prepared\n as a text file that contains two columns.\n
      1. Calculation of amplitude and exponent can\n be done for a data in eather lin-lin or log-log scales\n
      2. To get the amplitudes for a scaling form of\n the size characteristic data has to be in lin-lin scale"""
      tk.Label(newWindow,text=message).pack()

root=tk.Tk()
root.title("Scaling and such")
root.geometry('400x800')
root.resizable(0, 0)

###Building a scroll bar on the canvas
main_frame=tk.Frame(root)
main_frame.pack(fill='both',expand=1)

canvas=tk.Canvas(main_frame)
canvas.pack(side='left', fill='both',expand=1)

scrolbar=ttk.Scrollbar(main_frame,orient="vertical",command=canvas.yview)
scrolbar.pack(side="right",fill="y")

canvas.configure(yscrollcommand=scrolbar.set)
canvas.bind("<Configure>",lambda e:canvas.configure(scrollregion=canvas.bbox("all")))

wind=tk.Frame(canvas)
canvas.create_window((0,0),window=wind,anchor="nw")

### variables
scales=("lin-lin","log-log")
scale=tk.StringVar(value="lin-lin")
order=tk.StringVar(value="2")
start=tk.StringVar(value="0")
accuracy_Rg=tk.StringVar(value="")
amplitude=tk.StringVar(value="")
exponent=tk.StringVar(value="")
message_nu=tk.StringVar(value="")
message_amp=tk.StringVar(value="")
message_FS=tk.StringVar(value="")
message_RG=tk.StringVar(value="")
message_Scaling=tk.StringVar(value="")
N=tk.StringVar(value="1")
step=tk.StringVar(value="1")
step_g=tk.StringVar(value="100")
sc=tk.StringVar(value="")
scale.set("lin-lin")
order.set("2")
minsize=tk.StringVar(value="5000")
accuracy=tk.StringVar(value="0.05")
message_diff=tk.StringVar(value="")
num_mol=tk.StringVar(value="1")
message_auto=tk.StringVar(value="")
coef=tk.StringVar(value="")
score=tk.StringVar(value="")

title=tk.Label(wind,
                text="Size Scaling and diffution in polymers",
                font="Calibri 14 bold",anchor=tk.CENTER).pack(fill=tk.X)
tk.Label(wind,text=("This is a pack of functions that works with scaling\n in polymers. Please chose the data for processing\n (.csv,.txt or .dat)"),
                      font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
tk.Button(wind, text="Instructions", command=data_preparation).pack()


### Scaling exponent and amplitude
tk.Label(wind,
                text=("It allows to:\n"+
                "1. Get a calculation of amplitude and exponent \n for a data that is expected to have a scaling behavour"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
Frame_nu=tk.Frame(wind)
Frame_nu.columnconfigure(3)
Frame_nu.rowconfigure(1)
tk.OptionMenu(Frame_nu,scale,*scales).grid(row=0,column=0)
ttk.Button(Frame_nu, text="Calculate", command=lambda:calculate("Scaling_exponent_and_amplitude",(open_file(),scale.get()),start,accuracy_Rg,amplitude,exponent,message_nu)).grid(row=0,column=1)
ttk.Button(Frame_nu, text="Show the fit", command=lambda:display_picture("Scaling_exponent_and_amplitude",start)).grid(row=0,column=2)
Frame_nu.pack()
### Amplitudes for a full scaling form
tk.Label(wind,
                text=("2. Get a calculation of amplitudes for a size\n characteristic using known values of both size exponent \n and correction to scaling exponent"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
Frame_par=tk.Frame(wind)
Frame_par.columnconfigure(2)
Frame_par.rowconfigure(2)
tk.Label(Frame_par,
                text=("Please provide 1 if data has size characteristic and not its square"),
                font="Calibri 10",anchor=tk.CENTER).grid(row=0,column=0,columnspan=2)
tk.Entry(Frame_par,textvariable=order).grid(row=1,column=0)
ttk.Button(Frame_par, text="Calculate", command=lambda:calculate("Scaling_form_RG_chain",(open_file(),order.get()),amplitude,exponent,message_amp)).grid(row=1,column=1)
Frame_par.pack()

### Prediction of Rg for a given N

tk.Label(wind,
                text=("3. Get a prediction of gyration radius value \n for a linear chain with N monomers"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
Frame_pred=tk.Frame(wind)
Frame_pred.columnconfigure(3)
Frame_pred.rowconfigure(2)
Sim_type=tk.StringVar(value=None)
methods=("MD(build-in)","MC(build-in)","From my data")
Sim_type.set("From my data")
tk.OptionMenu(Frame_pred,Sim_type,*methods).grid(row=1,column=0)
tk.Label(Frame_pred,text="Number of monomers").grid(row=0,column=1)
tk.Entry(Frame_pred,textvariable=N).grid(row=1,column=1)
st=Sim_type.get()
def file(st):
    if st == "From my data": return open_file()
ttk.Button(Frame_pred, text="Calculate", command=lambda:calculate("Rg_chain_prediction",(N.get(),Sim_type.get(),file(Sim_type.get())),amplitude,exponent,message_RG)).grid(row=1,column=2)
Frame_pred.pack()

###finit size scaling
Frame_fin=tk.Frame(wind)
Frame_fin.columnconfigure(2)
Frame_fin.rowconfigure(1)
tk.Label(wind,
                text=("4. Finit size scaling calculation for\n size ratios and shape characteristics"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
ttk.Button(Frame_fin, text="Calculate", command=lambda:calculate("finit_size_scaling",(open_file()),start,amplitude,exponent,message_FS)).grid(row=0,column=0)
ttk.Button(Frame_fin, text="Show the fit", command=lambda:display_picture("finit_size_scaling",start)).grid(row=0,column=1)
Frame_fin.pack()
###Is there scaling
tk.Label(wind,
                text=("5. Checks if there is power law or\n exponential dependance in the tail of data"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
tk.Label(wind,text="A number of points added to the cutoff on each iteration step").pack()
tk.Entry(wind,textvariable=step).pack()
ttk.Button(wind, text="Calculate", command=lambda:calculate("Is_there_scaling",(open_file(),step.get()),amplitude,exponent,message_Scaling)).pack()


### diffution search
tk.Label(wind,
                text=("6. Checks if there is a diffusive behavior\n on long times for dynamic functions in melt"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)

tk.Label(wind,text="A number of points added to the cutoff on each iteration step").pack()
tk.Entry(wind,textvariable=step_g).pack()
tk.Label(wind,text="Minimum number of data points over which the fit is to be calculated").pack()
tk.Entry(wind,textvariable=minsize).pack()
tk.Label(wind,text="How close to 1 the exponent has to be").pack()
tk.Entry(wind,textvariable=accuracy).pack()
Frame_diff=tk.Frame(wind)
Frame_diff.columnconfigure(2)
Frame_diff.rowconfigure(1)
ttk.Button(Frame_diff, text="Calculate", command=lambda:calculate("g1_g3",(open_file(),step_g.get(),minsize.get(),accuracy.get()),start,coef,score,message_diff)).grid(row=0,column=0)
ttk.Button(Frame_diff, text="Show the fit", command=lambda:display_picture("g1_g3",start)).grid(row=0,column=1)
Frame_diff.pack()
'''
tk.Label(wind,
                text=("7. Calculates a relaxation time\n from autocorrelation function"),
                font="Calibri 12",anchor=tk.CENTER).pack(fill=tk.X)
tk.Label(wind,text="Number of molecules in the simulation box").pack()
tk.Entry(wind,textvariable=num_mol).pack()
ttk.Button(wind, text="Calculate", command=lambda:calculate("Autocorrelation",(open_file(),num_mol.get()),amplitude,exponent,message_auto)).pack()
tk.Label(wind,text="").pack()'''
root.mainloop()
        

