from setuptools import setup

setup(
    name="deblur_testing",
    version='0.1',
    py_modules=['script'],
    install_requires=[
        'qiime2 == 2017.11.*',
    ],
    entry_points='''
        [console_scripts]
        rdemux=script:demux
        rdeblur=script:pre_trims
        rpost_trim=script:post_trims
        ranalysis=script:analysis
    ''',
)

# command to install
# pip install --prefix /anaconda3/envs/qiime2-2017.11/ -e .