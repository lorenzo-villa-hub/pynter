import paramiko
import os

class PersistentSSH:
    def __init__(self, host, username=None, password=None, key_filename=None, port=None, timeout=10):
        """
        Initializes an SSH connection, using settings from ~/.ssh/config if available.

        :param host: The SSH host (can be an alias from ~/.ssh/config)
        :param username: SSH username (overrides config if provided)
        :param password: SSH password (optional)
        :param key_filename: Path to private key file (overrides config if provided)
        :param port: SSH port (overrides config if provided)
        :param timeout: Timeout for connection (default: 10 seconds)
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
        """Establishes an SSH connection using the loaded settings."""
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

    def execute_command(self, command, timeout=5):
        """Executes a command over SSH with a timeout."""
        if self.client:
            try:
                transport = self.client.get_transport()
                if transport and transport.is_active():
                    session = transport.open_session()
                    session.settimeout(timeout)
                    session.exec_command(command)
                    
                    stdout = session.makefile('r', -1).read()
                    stderr = session.makefile_stderr('r', -1).read()

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

# Example usage
if __name__ == "__main__":
    ssh = PersistentSSH("my-server-alias")  # Uses settings from ~/.ssh/config
    
    stdout, stderr = ssh.execute_command("ls -la")
    print("STDOUT:", stdout)
    print("STDERR:", stderr)

    ssh.close()
