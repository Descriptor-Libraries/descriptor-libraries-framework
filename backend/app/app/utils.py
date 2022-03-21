import subprocess


def obabel(in_format, out_format, data):
    proc = subprocess.run(f'echo "{data}" | obabel -ixyz -osdf', shell=True, capture_output=True)
    return proc.stdout.decode()
