
import os
import tempfile
import unittest 
from unittest.mock import patch, MagicMock, call
from pynter.hpc.ssh import SSHProtocol , get_ssh_config
from pynter.hpc.ssh import get_path_relative_to_hpc, get_rsync_to_hpc_command, get_rsync_from_hpc_command


class TestSSHConfig(unittest.TestCase):
    
    def test_get_ssh_config(self):
        # Create a temporary SSH config file
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            temp_file.write(
                "Host dummyhost\n"
                "    HostName example.com\n"
                "    Port 2222\n"
                "    User testuser\n"
                "    IdentityFile ~/.ssh/id_test\n"
            )
            temp_file_path = temp_file.name

        try:
            # Call the function with the path to the temp config
            host, port, username, key_filename = get_ssh_config(host="dummyhost",
                                                                port=22, #default value
                                                                path=temp_file_path)
            self.assertEqual(host, "example.com")
            self.assertEqual(port, 2222)
            self.assertEqual(username, "testuser")
            self.assertTrue(key_filename.endswith("id_test"))  # path expands later, so we just check the ending

        finally:
            os.remove(temp_file_path)  # Cleanup
 

def setup_mock_sshclient(MockSSHClient):
    mock_client = MagicMock()
    MockSSHClient.return_value = mock_client

    # Simulate active transport
    mock_transport = MagicMock()
    mock_transport.is_active.return_value = True
    mock_client.get_transport.return_value = mock_transport

    # Mock session and output
    mock_session = MagicMock()
    mock_transport.open_session.return_value = mock_session

    # stdout
    mock_stdout = MagicMock()
    mock_stdout.read.return_value = b"mocked test"
    mock_session.makefile.return_value = mock_stdout

    # stderr
    mock_stderr = MagicMock()
    mock_stderr.read.return_value = b""
    mock_session.makefile_stderr.return_value = mock_stderr

    return mock_client, mock_session


class TestSSHProtocol(unittest.TestCase):

    def setUp(self):
        patcher = patch("pynter.hpc.ssh.paramiko.SSHClient")
        self.addCleanup(patcher.stop)
        self.MockSSHClient = patcher.start()
        self.mock_client, self.mock_session = setup_mock_sshclient(self.MockSSHClient)    
   
    def test_args(self):
        ssh = SSHProtocol(host="fakehost",
                          username="user",
                          password='password',
                          key_filename='key_filename',
                          port=8820,
                          timeout=30,
                          reinitialize=True)
        
        self.assertEqual(ssh.host,'fakehost')
        self.assertEqual(ssh.username,'user')
        self.assertEqual(ssh.password,'password')
        self.assertEqual(ssh.key_filename,'key_filename')
        self.assertEqual(ssh.port,8820)
        self.assertEqual(ssh.timeout,30)
        self.assertEqual(ssh.initialized,True)

        self.mock_client.connect.assert_called_once_with(
              'fakehost',
              port=8820,
              username='user',
              password='password',
              key_filename='key_filename',
              timeout=30)

    def test_connect(self):
        ssh = SSHProtocol(host="fakehost", username="user", reinitialize=True)
        self.mock_client.connect.assert_called_once_with(
                    "fakehost",
                    port=22,
                    username="user",
                    password=None,
                    key_filename=None,
                    timeout=10)


    def test_run_command(self):
        ssh = SSHProtocol(host="fakehost", username="user", reinitialize=True)
        stdout, stderr = ssh.run_command("echo mocked test", printout=False)
        self.assertEqual(stdout, "mocked test")
        self.assertEqual(stderr, "")
        
        
    def test_instances_callback(self):
        ssh1 = SSHProtocol(host="fakehost", username="user", reinitialize=True)
        ssh2 = SSHProtocol(host="fakehost", username="user", reinitialize=False)
        ssh3 = SSHProtocol(host="fakehost", username="user", reinitialize=True,password='password')
        ssh4 = SSHProtocol(host="fakehost", username="new-user", reinitialize=False)
        ssh5 = SSHProtocol(host="new-fakehost", username="new-user", reinitialize=False)

        self.mock_client.connect.assert_has_calls([
            call("fakehost", port=22, username="user", password=None, key_filename=None, timeout=10),
            call("fakehost", port=22, username="user", password='password', key_filename=None, timeout=10),
            call("fakehost", port=22, username="new-user", password=None, key_filename=None, timeout=10),
            call("new-fakehost", port=22, username="new-user", password=None, key_filename=None, timeout=10)
            ],any_order=False)  

        instances = SSHProtocol._instances
        self.assertIn('fakehost',instances.keys())
        self.assertIn('user',instances['fakehost'].keys())
        self.assertIn('new-user',instances['fakehost'].keys())
        self.assertIn('new-fakehost',instances.keys())


class TestRsync(unittest.TestCase):
    
    def test_get_path_relative_to_hpc(self):
        # local path
        path = '/home/local/data/test/directory'
        localdir = '/home/local/data'
        remotedir = '/work/remote/data'
        path_in_hpc, path_relative, is_path_on_hpc = get_path_relative_to_hpc(path=path,
                                                                              localdir=localdir,
                                                                              remotedir=remotedir)
        self.assertEqual(path_in_hpc,'/work/remote/data/test/directory')
        self.assertEqual(path_relative,'/test/directory')
        self.assertEqual(is_path_on_hpc,False)

        # remote path
        path = '/work/remote/data/test/directory'
        localdir = '/home/local/data'
        remotedir = '/work/remote/data'
        path_in_hpc, path_relative, is_path_on_hpc = get_path_relative_to_hpc(path=path,
                                                                              localdir=localdir,
                                                                              remotedir=remotedir)
        self.assertEqual(path_in_hpc,path)
        self.assertEqual(path_relative,'/test/directory')
        self.assertEqual(is_path_on_hpc,True)   
        
        # path outside of localdir - should raise warning
        path = '/tmp/test/data/test/directory'
        localdir = '/home/local/data'
        remotedir = '/work/remote/data'
        path_in_hpc, path_relative, is_path_on_hpc = get_path_relative_to_hpc(path=path,
                                                                              localdir=localdir,
                                                                              remotedir=remotedir)
        self.assertEqual(is_path_on_hpc,False)   
        
        
    def test_rsync_from_hpc_command(self):
        
        def check_strings_in_command(strings,command):
            for string in strings:
                self.assertIn(string,command)
        
        command = get_rsync_from_hpc_command(hostname='fakehost',
                                             remotedir='/work/remote',
                                             localdir='/home/local',
                                             exclude=['file_ex1','file_ex2'],
                                             dry_run=True,
                                             makedir_on_local=True,
                                             rsync_args={
                                                 'files-from':'filelist.txt',
                                                 'relative':None
                                                 })
        desired_strings = ['mkdir -p /home/local/','--exclude=file_ex1', '--exclude=file_ex2',
                   '--files-from=filelist.txt','--relative','-e ssh fakehost:/work/remote/ /home/local/']
        check_strings_in_command(desired_strings, command)
       
        
        from pynter import LOCAL_DIR, REMOTE_DIR
        import os.path as op
        localdir = op.path.join(LOCAL_DIR,'data/directory')
        command = get_rsync_from_hpc_command(hostname='fakehost',localdir=localdir)        
        
        
        
        
        