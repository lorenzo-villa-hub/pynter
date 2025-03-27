
import os
import os.path as op
from shutil import which

from monty.dev import requires
import paramiko
from pynter import run_local, SETTINGS

class SSHProtocol:
    
    def __init__(self, host, username=None, password=None, key_filename=None, port=None, timeout=10):
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
        """
        self.host = host
        self.username = username
        self.password = password
        self.key_filename = key_filename
        self.port = port
        self.timeout = timeout
        self.client = None

        # Load SSH config if available
        self.load_ssh_config()

        # Connect using the loaded settings
        self.connect()


    def __enter__(self):
        """Enter the context: return the instance itself."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context: ensure the connection is closed."""
        self.close()


    def load_ssh_config(self):
        """Loads settings from ~/.ssh/config for the given host."""
        ssh_config_path = os.path.expanduser("~/.ssh/config")

        if os.path.exists(ssh_config_path):
            ssh_config = paramiko.config.SSHConfig()
            with open(ssh_config_path) as f:
                ssh_config.parse(f)

            # Get the host-specific settings (if defined)
            host_config = ssh_config.lookup(self.host)

            self.host = host_config.get("hostname", self.host)
            self.port = int(host_config.get("port", self.port or 22))
            self.username = self.username or host_config.get("user")
            self.key_filename = self.key_filename or host_config.get("identityfile", [None])[0]


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


    def run_command(self, command, timeout=10, printout=True):
        """Executes a command over SSH with a timeout."""
        if self.client:
            try:
                transport = self.client.get_transport()
                if transport and transport.is_active():
                    session = transport.open_session()
                    session.settimeout(timeout)
                    session.exec_command(command)
                    
                    stdout = session.makefile('r', -1).read().decode().strip() # decode string from b"" format
                    stderr = session.makefile_stderr('r', -1).read().decode().strip()
                    
                    if printout:
                        print("%s@%s:" %(self.username,self.host),command)
                        if stdout:
                            print(stdout)
                        if stderr:
                            print("ERROR:",stderr)

                    return stdout, stderr

            except Exception as e:
                return None, f"Command execution failed: {e}"
        else:
            return None, "SSH client is not connected"


    def close(self):
        """Closes the SSH connection."""
        if self.client:
            self.client.close()
            self.client = None
            


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
        Kwargs for the run_local function.

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
    run_local(localcmd)
    
    cmd = "rsync -r -uavzh " #keep the spaces
    if dry_run:
        cmd += "--dry-run "
    if exclude:
        for s in exclude:
            cmd += f'--exclude={s} ' 
    cmd += f"-e ssh {hostname}:{remotedir} {localdir} "

    stdout,stderr = run_local(cmd,printout=printout,dry_run=False,**kwargs)
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
        Kwargs for the run_local function.

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
    
    self.mkdir(remotedir,printout=False)
    
    cmd = "rsync -r -uavzh " #keep the spaces
    if dry_run:
        cmd += "--dry-run "
    if exclude:
        for s in exclude:
            cmd += f'--exclude={s} '
    cmd += f"-e ssh  {localdir} {hostname}:{remotedir} "
    
    stdout,stderr = run_local(cmd,printout=printout,dry_run=False,**kwargs)

    return stdout,stderr

