import paramiko

def test_server_connection(host, port, user, key_file):
    """
    Attempts to connect via SSH to the server and run a simple command.
    Prints the output if successful, or an error message if it fails.
    """
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        print(f"Connecting to {host}:{port} as {user} using key {key_file}")
        client.connect(hostname=host, port=port, username=user, key_filename=key_file)
        print("Connection successful!")

        # Run a test command
        stdin, stdout, stderr = client.exec_command("echo 'Hello from the server!'")
        output = stdout.read().decode().strip()
        err_output = stderr.read().decode().strip()

        print(f"Server response: {output}")
        if err_output:
            print(f"Server error (if any): {err_output}")

    except Exception as e:
        print(f"Connection failed: {str(e)}")
    finally:
        client.close()

if __name__ == "__main__":
    server = "203.101.229.234"
    port = 22
    user = "mfreeman"
    key_file = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"

    test_server_connection(server, port, user, key_file)
