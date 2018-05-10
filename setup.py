from setuptools import setup

setup(
    name="deblur_testing",
    version='0.1',
    py_modules=['scripts'],
    install_requires=[
        'Click',
        'qiime2 >= 2018.2.*',
        'wget'
    ],
    entry_points='''
        [console_scripts]
        rdemux=scripts:do_demux
        rdeblur=scripts:pre_trims
        rpost_trim=scripts:post_trims
        ranalysis=scripts:analysis
        rplot=scripts:do_plots
        rall=scripts:pre_post
        sbiom=scripts:biom_to_post
        biom2qza=scripts:biom_to_qiime
        qza2biom=scripts:qiime_to_biom
    ''',
)

# command to install
# pip install --prefix /anaconda3/envs/qiime2-2017.11/ -e .