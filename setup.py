from distutils.core import setup
setup(
    name='rbsa',
    version='v0.1',
    py_modules=['module.rbsa','module.ovlp2graph','module.graph2chr','module.chr_paths','module.filter'],
    author = 'likui',
    author_email = 'likui345@126.com',
    requires = ['mummer','networkx'],
    url='https://github.com/likui345/PGA'
)
