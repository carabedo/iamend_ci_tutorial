# iamend_ci

Libreria para la estimacion de la permeabilidad relativa efectiva usando mediciones de impedancia de una bobina sobre materiales conductores ferromagneticos. Contiene submodulos para la importancion de mediciones realizadas con `Solartron 1260A`, grafico de impedancias y correciones de efectos no ideales en los contactos.


### carga  y correccion datos



```python

data=cn.so.load('C:/Users/fernando/tesis/labs/impedancias/pp1/enero 18')
f=cn.so.getf(data)
bo=cn.bo.bobpp1
datacorr=cn.so.corr(f,bo,data)
```

