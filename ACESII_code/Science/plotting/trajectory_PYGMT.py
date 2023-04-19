from ACESII_code.class_var_func import setupPYGMT
setupPYGMT()
import pygmt
pygmt.show_versions()

fig = pygmt.Figure()
fig.coast(region="g", frame=True, shorelines=1)
fig.show()