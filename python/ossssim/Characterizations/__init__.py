import os

surveys = {}
for _ in os.listdir(__path__[0]):
    if _.startswith('_'):
        continue
    surveys[_]=os.path.join(__path__[0], _)

