from setuptools import setup

setup(
    name="deblur_testing",
    version='0.1',
    py_modules=['scripts'],
    install_requires=[
        'Click',
        'qiime2'
    ],
    entry_points='''
        [console_scripts]
        rdemux=scripts:do_demux
        rdeblur=scripts:pre_trims
        rpost_trim=scripts:post_trims
        ranalysis=scripts:analysis
        rplot=scripts:do_plots
        rall=scripts:pre_post
    ''',
)

# command to install
# pip install --prefix /anaconda3/envs/qiime2-2017.11/ -e .