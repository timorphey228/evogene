import unittest
from unittest.mock import patch
from pathlib import Path
from evogene.tools import find_tool, check_tools

class ExecutableDiscoveryTests(unittest.TestCase):
    @patch('evogene.tools.os.access',return_value=True)
    @patch('evogene.tools.Path.is_file',return_value=True)
    @patch('evogene.tools.shutil.which',return_value=None)
    @patch('evogene.tools.Path.home',return_value=Path('/Users/macbook'))
    def test_finds_datasets_in_home_bin(self,*_):
        self.assertEqual(find_tool('datasets'), '/Users/macbook/bin/datasets')

    @patch('evogene.tools.os.access',return_value=True)
    @patch('evogene.tools.Path.is_file',return_value=True)
    @patch('evogene.tools.shutil.which',return_value=None)
    def test_finds_minimap2_without_path(self,*_):
        self.assertTrue(find_tool('minimap2').endswith('/minimap2'))

    @patch('evogene.tools.find_tool',return_value=None)
    def test_preflight_missing(self,*_):
        result=check_tools()
        self.assertFalse(result['datasets']['ok'])
        self.assertFalse(result['minimap2']['ok'])

    @patch('evogene.tools.subprocess.run')
    @patch('evogene.tools.find_tool')
    def test_preflight_checks_versions_without_network(self, find, run):
        from subprocess import CompletedProcess
        find.side_effect=lambda name:'/usr/local/bin/'+name
        run.side_effect=lambda command, **kwargs: CompletedProcess(
            args=command,returncode=0,stdout='version 1',stderr='')
        status=check_tools()
        self.assertTrue(status['datasets']['ok'])
        self.assertTrue(status['minimap2']['ok'])
        self.assertEqual(
            [x.args[0] for x in run.call_args_list],
            [['/usr/local/bin/datasets','--version'],
             ['/usr/local/bin/minimap2','--version']]
        )

    def test_rejects_untrusted_name(self):
        with self.assertRaises(ValueError):
            find_tool('../../bad')
