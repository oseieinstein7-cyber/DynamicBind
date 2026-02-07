#!/usr/bin/env python3
"""
Polaris Connection Script for Windows

Connects to ALCF Polaris via SSH and allows file transfer.

Usage:
    python connect_polaris.py              # Interactive shell
    python connect_polaris.py --cmd "ls"   # Run single command
    python connect_polaris.py --upload local_file remote_path
    python connect_polaris.py --download remote_path local_file
"""

import paramiko
import getpass
import argparse
import os
from scp import SCPClient

# ALCF Polaris configuration
HOSTNAME = "polaris.alcf.anl.gov"
USERNAME = "eistien"  # Change to your username


def create_ssh_client(hostname, username, password):
    """Create and return an SSH client connection."""
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname, username=username, password=password)
    return ssh


def run_command(ssh, command):
    """Execute a command and return output."""
    stdin, stdout, stderr = ssh.exec_command(command)
    output = stdout.read().decode('utf-8')
    errors = stderr.read().decode('utf-8')
    return output, errors


def upload_file(ssh, local_path, remote_path):
    """Upload a file to Polaris."""
    with SCPClient(ssh.get_transport()) as scp:
        scp.put(local_path, remote_path)
    print(f"Uploaded: {local_path} -> {remote_path}")


def download_file(ssh, remote_path, local_path):
    """Download a file from Polaris."""
    with SCPClient(ssh.get_transport()) as scp:
        scp.get(remote_path, local_path)
    print(f"Downloaded: {remote_path} -> {local_path}")


def interactive_shell(ssh):
    """Start an interactive shell session."""
    channel = ssh.invoke_shell()

    print("=" * 60)
    print("Connected to Polaris! Type 'exit' to quit.")
    print("=" * 60)

    import sys
    import select
    import socket

    try:
        while True:
            # Get user input
            cmd = input("polaris> ")
            if cmd.lower() in ['exit', 'quit']:
                break

            # Send command
            channel.send(cmd + '\n')

            # Wait for response
            import time
            time.sleep(0.5)

            # Read output
            while channel.recv_ready():
                output = channel.recv(4096).decode('utf-8')
                print(output, end='')
    except KeyboardInterrupt:
        print("\nDisconnected.")


def main():
    parser = argparse.ArgumentParser(description="Connect to ALCF Polaris")
    parser.add_argument("--user", default=USERNAME, help="Username")
    parser.add_argument("--host", default=HOSTNAME, help="Hostname")
    parser.add_argument("--cmd", help="Run single command")
    parser.add_argument("--upload", nargs=2, metavar=("LOCAL", "REMOTE"),
                        help="Upload file")
    parser.add_argument("--download", nargs=2, metavar=("REMOTE", "LOCAL"),
                        help="Download file")

    args = parser.parse_args()

    # Get password
    print(f"Connecting to {args.user}@{args.host}")
    password = getpass.getpass("Password: ")

    try:
        ssh = create_ssh_client(args.host, args.user, password)
        print("Connected successfully!")

        if args.cmd:
            output, errors = run_command(ssh, args.cmd)
            if output:
                print(output)
            if errors:
                print("Errors:", errors)

        elif args.upload:
            upload_file(ssh, args.upload[0], args.upload[1])

        elif args.download:
            download_file(ssh, args.download[0], args.download[1])

        else:
            # Interactive mode
            interactive_shell(ssh)

        ssh.close()

    except paramiko.AuthenticationException:
        print("Authentication failed. Check username/password.")
    except paramiko.SSHException as e:
        print(f"SSH error: {e}")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
