## INTERACTIVE OUTLIERS
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox

def Fit_Iterativo(f,x,y,init,derr_1):
    pars=init
    i=0
    dist=len(init)
    while(dist>0.01*len(init) and i<10):
        pars0,covm=curve_fit(f,x,y,pars,derr_1)
        dist=(abs(pars-pars0)).sum()
        pars=pars0
        i=i+1
    return pars,covm

def SortbyX(x,dx,y,dy):
    a=x.argsort()
    X=np.array([])
    DX=X
    Y=X
    DY=X
    for i in a:
        X=np.append(X,x[i])
        DX=np.append(DX,dx[i])
        Y=np.append(Y,y[i])
        DY=np.append(DY,dy[i])
    return X,DX,Y,DY

class f0:
    def f(x,q):
        return q
    def der(x,q):
        return 0
    npars=1
    text="y = p0"

class f1:
    def f(x,m,q):
        return m*x+q
    def der(x,m,q):
        return m
    npars=2
    text="y = p0*x + p1"

class f2:
    def f(x,a,b,c):
        return a*x**2+b*x+c
    def der(x,a,b,c):
        return 2*a*x+b
    npars=3
    text="y = p0*x**2 + p1*x + p2"

class fE:
    def f(x,a,b,c):
        return a*np.exp(b*x)+c
    def der(x,a,b,c):
        return a*b*np.exp(b*x)
    npars=3
    text="y = p0*np.exp(p1*x) + p2"

fit_switch=f1
init=[0.,1.]
bars=0
err_eff=True

fig=plt.figure(1)
ax1=plt.subplot2grid((3,4), (0,0), rowspan=2, colspan=3)
ax1.grid(which='major')
ax1.minorticks_on()

ax2=plt.subplot2grid((3,4), (2,0), colspan=3, sharex=ax1)
ax2.grid(which='major')
ax2.minorticks_on()

x, dx, y, dy = SortbyX(*np.loadtxt('test.txt',unpack=True))

init,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,dy)
derr=np.array([])
if(err_eff):
    derr=np.sqrt(dy**2+(dx*fit_switch.der(x,*init))**2)
else:
    derr=dy
visible=[True for i in x]

pl1,aaa,bbb=ax1.errorbar(x,y,dy,dx,linestyle='',marker='.',color='red')
pl2,ccc,ddd=ax1.errorbar(x,y,dy,dx,linestyle='',marker='.',color='green')

ics=np.linspace(np.min(x)-0.1*(np.max(x)-np.min(x)),np.max(x)+0.1*(np.max(x)-np.min(x)),200)
pl0,=ax1.plot(ics,fit_switch.f(ics,*init),linewidth=1,color='blue')

pl3,=ax2.plot(x,(y-fit_switch.f(x,*init))/derr,linewidth=1,linestyle='dashed',color='purple',marker='.')
pl4,=ax2.plot(x,(y-fit_switch.f(x,*init))/derr,linewidth=1,linestyle='dashed',color='blue',marker='.')

def Ripeti(x,dx,y,dy,init):
    if(err_eff):
        derr_1=np.sqrt(dy**2+(dx*fit_switch.der(x,*init))**2)
    else:
        derr_1=dy
    popt,pcov=Fit_Iterativo(fit_switch.f,x,y,init,derr_1)
    pl2.set_data(x,y)
    errorx=np.array([np.array([[x[i]-dx[i],y[i]],[x[i]+dx[i],y[i]]]) for i in range(len(x))])
    errory=np.array([np.array([[x[i],y[i]-dy[i]],[x[i],y[i]+dy[i]]]) for i in range(len(x))])
    ddd[0].set_segments(errorx)
    ddd[1].set_segments(errory)
    pl0.set_data(ics,fit_switch.f(ics,*popt))
    pl4.set_data(x,(y-fit_switch.f(x,*popt))/derr_1)
    ax2.relim()
    ax2.autoscale_view()
    fig.canvas.flush_events()
    fig.canvas.draw()
    return popt,pcov

dragx=0.
dragy=0.
cliccato=False

ax3=plt.subplot2grid((4,4), (0,3), rowspan=2)
ax3.set_axis_off()
result_box=ax3.text(0,1,"Result:"+'\n',horizontalalignment='left',verticalalignment='top',multialignment='left',wrap=True,fontsize=8)

def get_index(xd,yd,xc,yc):
    mind=0.0004
    mini=-1
    for i in range(len(xd)):
        if(((xd[i]-xc)/np.max(xd))**2+((yd[i]-yc)/np.max(yd))**2<mind):
            mind=((xd[i]-xc)/np.max(xd))**2+((yd[i]-yc)/np.max(yd))**2
            mini=i
    return mini

def get_visible():
    dist=(fit_switch.f(x,*init)-y)/derr
    a=np.array([])
    b=np.array([])
    c=np.array([])
    d=np.array([])
    for i in range(len(x)):
        if(visible[i] and ((bars>0 and dist[i]**2<=bars**2) or bars<=0)):
            a=np.append(a,x[i])
            b=np.append(b,dx[i])
            c=np.append(c,y[i])
            d=np.append(d,dy[i])
    return a,b,c,d

def result():
    s=fit_switch.text+'\n'+'\n'
    s+="all data:"+'\n'
    for i in range(fit_switch.npars):
        s+="p%d = %f +- %f    " % (i,popt_0[i],pcov_0.diagonal()[i])
    s+='\n'
    if(err_eff):
        derr=np.sqrt(dy**2+(dx*fit_switch.der(x,*popt_0))**2)
    else:
        derr=dy
    chisq_0=((y-fit_switch.f(x,*popt_0))**2/derr**2).sum()
    s+="chisq = %f" % (chisq_0)+'\n'+'\n'
    X,dX,Y,dY=get_visible()
    s+="selected data"
    if(bars!=0):
        s+=" (within %d error bars)" %(bars)
    s+=":"+'\n'
    for i in range(fit_switch.npars):
        s+="p%d = %f +- %f    " % (i,popt_1[i],pcov_1.diagonal()[i])
    s+='\n'
    if(err_eff):
        dERR=np.sqrt(dY**2+(dX*fit_switch.der(X,*popt_1))**2)
    else:
        dERR=dy
    chisq_1=((Y-fit_switch.f(X,*popt_1))**2/dERR**2).sum()
    s+="chisq = %f" % (chisq_1)+'\n'
    return s

def up_drag(event):
    global dragx,dragy
    dragx=event.xdata
    dragy=event.ydata
    cliccato=True

def on_click(event):
    global dragx, dragy, popt_1, pcov_1
    if(event.key=='shift'):
        reset()
        return
    if(not event.xdata or not event.ydata):
        return
    if((dragx-event.xdata)**2>0.01*ax1.get_xlim()[1] or (dragy-event.ydata)**2>0.01*ax1.get_ylim()[1]):
        return
    if(event.button==1 or event.button==3):
        i=get_index(x,y,event.xdata,event.ydata)
        if(i==-1):
            return
        if(visible[i]):
            visible[i]=False
        else:
            visible[i]=True
        X,dX,Y,dY=get_visible()
        popt_1,pcov_1=Ripeti(X,dX,Y,dY,init)
        result_box.set_text("Result:"+'\n'+result())
        fig.canvas.flush_events()
        fig.canvas.draw()
    return

def on_click_plus(event):
    on_click(event)
    cliccato=False
    return

def on_drag(event):
    global dragx,dragy
    if(not(cliccato)):
        return
    if(((dragx-event.xdata)**2<0.01*ax1.get_xlim()[1] and (dragy-event.ydata)**2<0.01*ax1.get_ylim()[1])):
        return
    dragx=event.xdata
    dragy=event.ydata
    on_click(event)
    return

def reset():
    global visible, dragx, dragy, bars, popt_1, pcov_1
    visible=[True for i in x]
    popt_1,pcov_1=Ripeti(x,dx,y,dy,init)
    result_box.set_text("Result:"+'\n'+result())
    dragx=0.
    dragy=0.
    bars=0
    return

def on_key(event):
    global fit_switch, init, bars, popt_0, pcov_0, popt_1, pcov_1, derr

    if(event.key=='='):
        fit_switch=f0
        init=[0.]
        init,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,dy)
        if(err_eff):
            derr=np.sqrt(dy**2+(dx*fit_switch.der(x,*init))**2)
        else:
            derr=dy
        popt_0,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,derr)
        init=popt_0
        pl5.set_data(ics,fit_switch.f(ics,*init))
        pl3.set_data(x,(y-fit_switch.f(x,*init))/derr)
        X,dX,Y,dY=get_visible()
        popt_1,pcov_1=Ripeti(X,dX,Y,dY,init)
        result_box.set_text("Result:"+'\n'+result())
        fig.canvas.flush_events()
        fig.canvas.draw()

    elif(event.key=='!'):
        fit_switch=f1
        init=[0.,1.]
        init,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,dy)
        if(err_eff):
            derr=np.sqrt(dy**2+(dx*fit_switch.der(x,*init))**2)
        else:
            derr=dy
        popt_0,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,derr)
        init=popt_0
        pl5.set_data(ics,fit_switch.f(ics,*init))
        pl3.set_data(x,(y-fit_switch.f(x,*init))/derr)
        X,dX,Y,dY=get_visible()
        popt_1,pcov_1=Ripeti(X,dX,Y,dY,init)
        result_box.set_text("Result:"+'\n'+result())
        fig.canvas.flush_events()
        fig.canvas.draw()

    elif(event.key=='"'):
        fit_switch=f2
        init=[0.,1.,2.]
        init,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,dy)
        if(err_eff):
            derr=np.sqrt(dy**2+(dx*fit_switch.der(x,*init))**2)
        else:
            derr=dy
        popt_0,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,derr)
        init=popt_0
        pl5.set_data(ics,fit_switch.f(ics,*init))
        pl3.set_data(x,(y-fit_switch.f(x,*init))/derr)
        X,dX,Y,dY=get_visible()
        popt_1,pcov_1=Ripeti(X,dX,Y,dY,init)
        result_box.set_text("Result:"+'\n'+result())
        fig.canvas.flush_events()
        fig.canvas.draw()

    elif(event.key=='E'):
        fit_switch=fE
        init=[0.,1.,2.]
        init,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,dy)
        if(err_eff):
            derr=np.sqrt(dy**2+(dx*fit_switch.der(x,*init))**2)
        else:
            derr=dy
        popt_0,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,derr)
        init=popt_0
        pl5.set_data(ics,fit_switch.f(ics,*init))
        pl3.set_data(x,(y-fit_switch.f(x,*init))/derr)
        X,dX,Y,dY=get_visible()
        popt_1,pcov_1=Ripeti(X,dX,Y,dY,init)
        result_box.set_text("Result:"+'\n'+result())
        fig.canvas.flush_events()
        fig.canvas.draw()

    elif(event.key=='1' or event.key=='2' or event.key=='3' or event.key=='4' or event.key=='5' or event.key=='0'):
        bars=eval(event.key)
        X,dX,Y,dY=get_visible()
        popt_1,pcov_1=Ripeti(X,dX,Y,dY,init)
        result_box.set_text("Result:"+'\n'+result())
        fig.canvas.flush_events()
        fig.canvas.draw()


connid = fig.canvas.mpl_connect('button_release_event', on_click_plus)
connid = fig.canvas.mpl_connect('button_press_event', up_drag)
connid = fig.canvas.mpl_connect('motion_notify_event', on_drag)
connid = fig.canvas.mpl_connect('key_press_event', on_key)

popt_1,pcov_1=Ripeti(x,dx,y,dy,init)
popt_0,pcov_0=Fit_Iterativo(fit_switch.f,x,y,init,derr)
init=popt_0
ics=np.linspace(np.min(x)-0.1*(np.max(x)-np.min(x)),np.max(x)+0.1*(np.max(x)-np.min(x)),200)
pl5,=ax1.plot(ics,fit_switch.f(ics,*popt_0),linewidth=1,color='purple')
pl3.set_data(x,(y-fit_switch.f(x,*popt_0))/derr)
ax2.relim()
ax2.autoscale_view()
result_box.set_text("Result:"+'\n'+result())

plt.show()