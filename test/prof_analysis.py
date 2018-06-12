from click.testing import CliRunner
from scripts import analysis
import os

runner = CliRunner()
wd = os.getcwd() + "/data/prof_analysis/"
with runner.isolated_filesystem():
    result = runner.invoke(analysis,        ["-i", wd,
                                             "-o", wd,
                                             "-tl", 100,
                                             "-tl", 90])

    print(result.exc_info)
    print(result.output)
