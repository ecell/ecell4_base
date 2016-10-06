
# World and Simulator with ODE solver

If you read through [Introduction](http://nbviewer.ipython.org/github/ecell/ecell4/blob/master/ipynb/Tutorials/Introduction.ipynb), it is NOT difficult to use **World** and **Simulator**.
**volume** and **{'C':60}** is equivalent of the **World** and **solver** is the **Simulator** below.


```python
%matplotlib inline
import numpy
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

y = run_simulation(numpy.linspace(0, 10, 100), {'C': 60}, volume=1.0)

```


![png](worldsim_1_0.png)


Here we give you a breakdown for **run_simulation**.
**run_simulation** use ODE simulator by default, so we create ODE world step by step.

## Creating ODE world

You can create world like this.


```python
w = ode.ODEWorld(Real3(1, 1, 1))
```

**Real3** is a coordinate vector.
In this example, the first argument for ODEWorld constructor is a cube.
Note that you can NOT use volume for ode.ODEWorld argument, like **run_simulation** argument.

Now you created a cube box for simulation, next let's throw molecules into the cube.


```python
w = ode.ODEWorld(Real3(1, 1, 1))
w.add_molecules(Species('C'), 60)
print(w.t(), w.num_molecules(Species('C')))  # will return (0.0, 60)
```

    (0.0, 60)


Use **add_molecules** to add molecules, **remove_molecules** to remove molecules, **num_molecules** to know the number of molecules.
First argument for each method is the **Species** you want to know.
You can get current time by **t** method.
However the number of molecules in ODE solver is real number, in these **_molecules** functions work only for integer number.
If you use real number in ODE, use **set_value** and **get_value**.


## How to use Real3

Before the detail of **Simulator**, we explaing more about **Real3**.


```python
pos = Real3(1, 2, 3)
print(pos)  # will print <ecell4.core.Real3 object at 0x7f44e118b9c0>
print(tuple(pos))  # will print (1.0, 2.0, 3.0)
```

    <ecell4.core.Real3 object at 0x7fc73ff9a948>
    (1.0, 2.0, 3.0)


You can not print **Real3** object directly.
You need to convert **Real3** to Python tuple or list once.


```python
pos1 = Real3(1, 1, 1)
x, y, z = pos[0], pos[1], pos[2]
pos2 = pos1 + pos1
pos3 = pos1 * 3
pos4 = pos1 / 5
print(length(pos1))  # will print 1.73205080757
print(dot_product(pos1, pos3))  # will print 9.0
```

    1.73205080757
    9.0


You can use basic function like dot_product.
Of course you can convert Real3 to numpy array too.


```python
a = numpy.asarray(tuple(Real3(1, 2, 3)))
print(a)  # will print [ 1.  2.  3.]
```

    [ 1.  2.  3.]


## Creating and Running ODESimulator

You can create a Simulator with Model and World like


```python
with reaction_rules():
    A + B > C | 0.01  # equivalent to create_binding_reaction_rule
    C > A + B | 0.3   # equivalent to create_unbinding_reaction_rule

m = get_model()

sim = ode.ODESimulator(m, w)
sim.run(10.0)
```

then call **run** method, the simulation will run.
In this example the simulation runs for 10seconds.

You can check the state of the **World** like this.


```python
print(w.t(), w.num_molecules(Species('C')))  # will return (10.0, 30)
```

    (10.0, 30)


You can see that the number of the **Species** **C** decreases from 60 to 30.

**World** describes the state of a timepoint, so you can NOT see the transition of the simulation with the **World**.
To obtain the time-series result, use **Observer**.


```python
obs = FixedIntervalNumberObserver(0.1, ('A', 'C'))
sim.run(10.0, obs)
print(obs.data())  # will return [[0.0, 0.0, 60.0], ..., [10.0, 29.994446899698026, 30.005553100301967]]
```

    [[10.0, 29.994445225953793, 30.005554774046203], [10.1, 29.994923291781845, 30.00507670821815], [10.2, 29.995360215494838, 30.00463978450516], [10.3, 29.99575953741398, 30.004240462586015], [10.4, 29.996124493269757, 30.00387550673024], [10.5, 29.996458040397858, 30.00354195960214], [10.6, 29.99676288168368, 30.003237118316317], [10.7, 29.99704148744875, 30.002958512551245], [10.8, 29.99729611545562, 30.002703884544378], [10.9, 29.997528829192905, 30.00247117080709], [11.0, 29.997741514588153, 30.002258485411843], [11.1, 29.997935895283558, 30.00206410471644], [11.2, 29.99811354659805, 30.001886453401948], [11.3, 29.998275908288566, 30.00172409171143], [11.4, 29.99842429621383, 30.001575703786166], [11.5, 29.998559912994793, 30.001440087005204], [11.6, 29.998683857758166, 30.00131614224183], [11.7, 29.99879713504173, 30.001202864958266], [11.8, 29.99890066293361, 30.001099337066385], [11.9, 29.998995280511267, 30.00100471948873], [12.0, 29.999081754640535, 30.00091824535946], [12.1, 29.999160786189673, 30.000839213810323], [12.2, 29.999233015708754, 30.000766984291243], [12.3, 29.999299028620428, 30.00070097137957], [12.4, 29.999359359963993, 30.000640640036004], [12.5, 29.999414498731316, 30.00058550126868], [12.6, 29.999464891829593, 30.000535108170403], [12.7, 29.999510947703186, 30.00048905229681], [12.8, 29.99955303964374, 30.000446960356257], [12.9, 29.999591508815463, 30.000408491184533], [13.0, 29.999626667020085, 30.000373332979912], [13.1, 29.999658799223873, 30.000341200776123], [13.2, 29.99968816586716, 30.000311834132837], [13.3, 29.999715004975137, 30.00028499502486], [13.4, 29.999739534087013, 30.000260465912984], [13.5, 29.999761952019092, 30.000238047980904], [13.6, 29.999782440476185, 30.00021755952381], [13.7, 29.9998011655243, 30.000198834475697], [13.8, 29.99981827893661, 30.000181721063388], [13.9, 29.999833919423562, 30.000166080576435], [14.0, 29.99984821375716, 30.000151786242835], [14.100000000000001, 29.99986127779844, 30.000138722201555], [14.2, 29.99987321743654, 30.000126782563456], [14.3, 29.99988412944696, 30.000115870553035], [14.4, 29.999894102275928, 30.00010589772407], [14.5, 29.999903216757303, 30.000096783242693], [14.600000000000001, 29.999911546767716, 30.00008845323228], [14.7, 29.999919159825385, 30.00008084017461], [14.8, 29.99992611763737, 30.000073882362628], [14.9, 29.99993247659969, 30.000067523400308], [15.0, 29.99993828825447, 30.000061711745527], [15.100000000000001, 29.999943599707688, 30.00005640029231], [15.2, 29.999948454010994, 30.000051545989002], [15.3, 29.999952890510645, 30.00004710948935], [15.4, 29.99995694516643, 30.000043054833565], [15.5, 29.999960650843132, 30.000039349156864], [15.600000000000001, 29.999964037576916, 30.00003596242308], [15.7, 29.999967132818774, 30.000032867181222], [15.8, 29.999969961657037, 30.00003003834296], [15.9, 29.999972547020704, 30.000027452979293], [16.0, 29.999974909865326, 30.00002509013467], [16.1, 29.999977069342826, 30.00002293065717], [16.2, 29.99997904295676, 30.000020957043237], [16.3, 29.99998084670416, 30.000019153295835], [16.4, 29.999982495205234, 30.000017504794762], [16.5, 29.999984001821836, 30.00001599817816], [16.6, 29.999985378765786, 30.00001462123421], [16.7, 29.99998663719784, 30.000013362802157], [16.8, 29.999987787318176, 30.00001221268182], [16.9, 29.99998883844905, 30.000011161550947], [17.0, 29.99998979911036, 30.000010200889637], [17.1, 29.999990677088707, 30.00000932291129], [17.2, 29.999991479500515, 30.00000852049948], [17.3, 29.99999221284971, 30.000007787150288], [17.4, 29.999992883080424, 30.000007116919573], [17.5, 29.999993495625187, 30.00000650437481], [17.6, 29.999994055448955, 30.00000594455104], [17.7, 29.99999456708936, 30.000005432910637], [17.8, 29.99999503469349, 30.000004965306506], [17.9, 29.999995462051498, 30.0000045379485], [18.0, 29.99999585262731, 30.000004147372685], [18.1, 29.99999620958673, 30.000003790413267], [18.200000000000003, 29.99999653582308, 30.000003464176917], [18.3, 29.999996833980653, 30.000003166019344], [18.4, 29.99999710647616, 30.000002893523835], [18.5, 29.99999735551831, 30.000002644481686], [18.6, 29.999997583125694, 30.000002416874302], [18.700000000000003, 29.999997791143187, 30.00000220885681], [18.8, 29.99999798125686, 30.000002018743135], [18.9, 29.99999815500768, 30.000001844992315], [19.0, 29.999998313803975, 30.000001686196022], [19.1, 29.999998458932858, 30.00000154106714], [19.200000000000003, 29.999998591570673, 30.000001408429323], [19.3, 29.99999871279251, 30.000001287207485], [19.4, 29.999998823580924, 30.000001176419072], [19.5, 29.999998924833914, 30.000001075166082], [19.6, 29.999999017372176, 30.00000098262782], [19.700000000000003, 29.999999101945786, 30.00000089805421], [19.8, 29.99999917924024, 30.000000820759755], [19.9, 29.999999249882055, 30.00000075011794], [20.0, 29.999999314443816, 30.00000068555618]]


There are several types of **Observer**s for E-Cell4.
**FixedIntervalNumberObserver** is the simplest **Observer** to obtain the time-series result.
As its name suggests, this **Observer** records the number of molecules for each time-step.
The 1st argument is the time-step, the 2nd argument is the molecule types.
You can check the result with **data** method, but there is a shortcut for this


```python
viz.plot_number_observer(obs)
```


![png](worldsim_24_0.png)


This plots the time-series result easily.

We explained the internal of **run_simulation** function.
When you change the **World** after creating the **Simulator**, you need to indicate it to **Simulator**.
So do NOT forget to call


```python
sim.initialize()
```

## Switching the solver
It is NOT difficult to switch the solver to stochastic method, as we showed **run_simulation**.


```python
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

m = get_model()

# ode.ODEWorld -> gillespie.GillespieWorld
w = gillespie.GillespieWorld(Real3(1, 1, 1))
w.add_molecules(Species('C'), 60)

# ode.ODESimulator -> gillespie.GillespieSimulator
sim = gillespie.GillespieSimulator(m, w)
obs = FixedIntervalNumberObserver(0.1, ('A', 'C'))
sim.run(10.0, obs)

viz.plot_number_observer(obs)
```


![png](worldsim_28_0.png)


**World** and **Simulator** do NOT change the **Model**, so you can switch several simulators for 1 model.
