
import os
import os.path as op
from shutil import which

from monty.dev import requires
import paramiko
from collections import defaultdict

from pynter import SETTINGS, run_command

def get_ssh_config(host,port=22,path=None):
    """
    Gets settings from config file for the given host.
    
    Parameters
    ----------
    host: (str) 
        The SSH host (can be an alias from ~/.ssh/config)
    port: (int)
        SSH port (default is 22).
    path: (str)
        Path of configuration file, if not provided uses the default ~/.ssh/config  
            
    Returns
    -------
    host: (str) 
        The SSH host.
    port: (int)
        SSH port.
    username: (str)
        SSH username (overrides config if provided)
    password: (str)
        SSH password.
    key_filename: (str)
        Path to private key file (overrides config if provided)
    """
    if path:
        path = os.path.abspath(path)
    ssh_config_path = path or os.path.expanduser("~/.ssh/config")

    if os.path.exists(ssh_config_path):
        ssh_config = paramiko.config.SSHConfig()
        with open(ssh_config_path) as f:
            ssh_config.parse(f)

        # Get the host-specific settings (if defined)
        host_config = ssh_config.lookup(host)

        host = host_config.get("hostname", host)
        port = int(host_config.get("port", port))
        username = host_config.get("user")
        key_filename = host_config.get("identityfile", [None])[0]

        return host, port, username, key_filename



class SSHProtocol:
    
    _instances = defaultdict(dict) #define attribute at class level, avoid having to define nested dict
    
    def __new__(cls, host, username=None, password=None, key_filename=None, port=22, timeout=10, reinitialize=False):
        """
        Creates or retrieves an SSHProtocol instance based on host and username.
        Ensures that when the same connection is not opened multiple times.
        If username is not provided, loads settings from  ~/.ssh/config.
        """
        if not username:
            host, port, username, key_filename = get_ssh_config(host)
        if username not in cls._instances[host] or reinitialize: 
            instance = super().__new__(cls)
            instance.host = host
            instance.username = username
            instance.password = password
            instance.key_filename = key_filename
            instance.port = port
            instance.timeout = timeout
            cls._instances[host][username] = instance
        return cls._instances[host][username]
    
    
    def __init__(self, 
                 host, 
                 username=None, 
                 password=None, 
                 key_filename=None, 
                 port=22, 
                 timeout=10, 
                 reinitialize=False):
        """
        Initializes an SSH connection, using settings from ~/.ssh/config if available.
        
        Parameters
        ----------

        host: (str) 
            The SSH host (can be an alias from ~/.ssh/config)
        username: (str)
            SSH username (overrides config if provided)
        password: (str)
            SSH password.
        key_filename: (str)
            Path to private key file (overrides config if provided)
        port: (int)
            SSH port (overrides config if provided)
        timeout: (int)
            Timeout for connection (default: 10 seconds)
        reinitialize: (bool)
            Force to create a new connection instead of using existing ones.
        """
        if not hasattr(self, 'initialized') or reinitialize:
            self.client = None
            self.initialized = True
            self.connect()
        else:
            if not self.client:
                print(f"Re-trying connection to {self.host}")
                self.connect()
                

    def __enter__(self):
        """Enter the context: return the instance itself."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context: ensure the connection is closed."""
        self.close()


    def __repr__(self):
        lines = [f"{key}: {getattr(self,key)}" for key in ['host','port','username']]
        return '\n'.join(lines)
    
    def __print__(self):
        return self.__repr__()


    def close(self):
        """Closes the SSH connection."""
        if self.client:
            self.client.close()
            self.client = None


    def connect(self):
        """
        Establishes an SSH connection using the loaded settings, reuses if already connected.
        """
        if self.client and self.client.get_transport() and self.client.get_transport().is_active():
            print(f"Reusing existing connection to {self.host}")
            return 
        self.client = paramiko.SSHClient()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        try:
            self.client.connect(
                self.host,
                port=self.port,
                username=self.username,
                password=self.password,
                key_filename=self.key_filename,
                timeout=self.timeout
            )
        except Exception as e:
            print(f"Failed to connect: {e}")
            self.client = None


    def load_ssh_config(self):
        """
        Loads settings from ~/.ssh/config for the given host.
        """
        host, port, username, key_filename = get_ssh_config(self.host,self.port)
        self.host = host
        self.port = port
        self.username = username
        self.key_filename = key_filename


    def run_command(self,
                    command, 
                    timeout=10, 
                    printout=True, 
                    dry_run=False):
        """
        Executes a command over SSH with a timeout.
        """
        message = "%s@%s:" %(self.username,self.host),command
        if dry_run:
            if printout:
                print(message)
            return command, ''
        elif self.client:
            try:
                transport = self.client.get_transport()
                if transport and transport.is_active():
                    session = transport.open_session()
                    session.settimeout(timeout)
                    session.exec_command(command)
                    
                    stdout = session.makefile('r', -1).read().decode().strip() # decode string from b"" format
                    stderr = session.makefile_stderr('r', -1).read().decode().strip()
                    
                    if printout:
                        print(message)
                        if stdout:
                            print(stdout)
                        if stderr:
                            print("ERROR:",stderr)

                    return stdout, stderr

            except Exception as e:
                return None, f"Command execution failed: {e}"
        else:
            return None, "SSH client is not connected"


            


@requires(which("rsync"),
      "rsync needs to be installed, you can install it with 'sudo apt-get install rsync'.")

def rsync_from_hpc(
                    hostname=None,
                    remotedir=None,
                    localdir=None,
                    exclude=None,
                    printout=True,
                    dry_run=False,
                    **kwargs
                    ):
    """
    Sync folders from HPC to local machine. The command "rsync" is used. With this function all
    the folders contained in the remote dir are synched to the local dir.

    Parameters
    ----------
    hostname : (str)
        Name of the host.
    remotedir : (str)
        Remote directory. The default is None. If None the one in config.yml file is used.
    localdir : (str)
        Local directory. The default is None. If None the one in config.yml file is used.
    exclude : (list)
        List of files to exclude in rsync. The default is None.
    dry_run : (bool)
        Perform dry run in rsync with --dry-run. The default is False. The dry_run in run_local is set to False.
    kwargs : (dict)
        Kwargs for the local run_command function.

    Returns
    -------
    stdout : (str)
        Output.
    stderr : (str)
        Error.
    """        
    hostname = hostname if hostname else SETTINGS['HPC']['hostname']
    remotedir = remotedir if remotedir else SETTINGS['HPC']['remotedir']
    localdir = localdir if localdir else SETTINGS['HPC']['localdir']
    localdir = op.abspath(localdir)
    remotedir = op.join(remotedir,'')  #ensure backslash at the end
    localdir = op.join(localdir,'')
    
    localcmd = 'mkdir -p %s' %localdir
    run_command(localcmd)
    
    cmd = "rsync -r -uavzh " #keep the spaces
    if dry_run:
        cmd += "--dry-run "
    if exclude:
        for s in exclude:
            cmd += f'--exclude={s} ' 
    cmd += f"-e ssh {hostname}:{remotedir} {localdir} "

    stdout,stderr = run_command(cmd,printout=printout,dry_run=False)
    return stdout,stderr


@requires(which("rsync"),
      "rsync needs to be installed, you can install it with 'sudo apt-get install rsync'.")

def rsync_to_hpc(
                hostname=None,
                localdir=None,
                remotedir=None,
                exclude=None,
                printout=True,
                dry_run=False,
                **kwargs
                ):
    """
    Sync folders from local machine to HPC. The command "rsync" is used. With this function all
    the folders contained in the local dir are synched to the remote dir.

    Parameters
    ----------
    hostname : (str)
        Name of the host.
    localdir : (str)
        Local directory. The default is None. If None the one in config.yml file is used.
    remotedir : (str)
        Remote directory. The default is None. If None the one in config.yml file is used.
    exclude : (list)
        List of files to exclude in rsync. The default is None.
    dry_run : (bool)
        Perform dry run in rsync with --dry-run. The default is False. The dry_run in run_local is set to False.
    kwargs : (dict)
        Kwargs for the local run_command function.

    Returns
    -------
    stdout : (str)
        Output.
    stderr : (str)
        Error.
    """        
    hostname = hostname if hostname else SETTINGS['HPC']['hostname']
    remotedir = remotedir if remotedir else SETTINGS['HPC']['remotedir']
    localdir = localdir if localdir else SETTINGS['HPC']['localdir']
    localdir = op.abspath(localdir)
    localdir = op.join(localdir,'')  #ensure backslash at the end
    remotedir = op.join(remotedir,'')
    
    ssh = SSHProtocol(hostname)
    ssh.run_command('mkdir -p %s' %remotedir)
    
    cmd = "rsync -r -uavzh " #keep the spaces
    if dry_run:
        cmd += "--dry-run "
    if exclude:
        for s in exclude:
            cmd += f'--exclude={s} '
    cmd += f"-e ssh  {localdir} {hostname}:{remotedir} "
    
    stdout,stderr = run_command(cmd,printout=printout,dry_run=False,**kwargs)

    return stdout,stderr

