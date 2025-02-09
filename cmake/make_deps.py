# SPDX-FileCopyrightText: 2024 Howetuft
#
# SPDX-License-Identifier: Apache-2.0

"""This script downloads and installs dependencies for LuxCore build."""

import os
import tempfile
from urllib.request import urlretrieve
from pathlib import Path
from zipfile import ZipFile
import subprocess
import logging
import shutil
import textwrap
import argparse
import json
import configparser
import platform

CONAN_APP = None

CONAN_ALL_PACKAGES = '"*"'

BUILD_DIR = os.getenv("BUILD_DIR", "./build")

CMAKE_DIR = Path(BUILD_DIR) / "cmake"

logger = logging.getLogger("LuxCoreDeps")

CONAN_ENV = {}

URL_SUFFIXES = {
    "Linux-X64": "ubuntu-latest",
    "Windows-X64": "windows-latest",
    "macOS-ARM64": "macos-14",
    "macOS-X64": "macos-13",
}

def find_platform():
    system = platform.system()
    if system == "Linux":
        res = "Linux-X64"
    elif system == "Windows":
        res = "Windows-X64"
    elif system == "Darwin":
        machine = platform.machine()
        if machine == "arm64":
            res = "macOS-ARM64"
        elif machine == "x86_64":
            res = "macOS-X64"
        else:
            raise RuntimeError(f"Unknown machine for MacOS: '{machine}'")
    else:
        raise RuntimeError(f"Unknown system '{system}'")
    return res


def build_url(user, release):
    """Build the url to download from."""
    suffix = URL_SUFFIXES[find_platform()]

    if not user:
        user = "LuxCoreRender"

    return f"https://github.com/{user}/LuxCoreDeps/releases/download/{release}/luxcore-deps-{suffix}.zip"


def get_profile_name():
    """Get the profile file name, based on platform."""
    return f"conan-profile-{find_platform()}"


def ensure_conan_app():
    logger.info("Looking for conan")
    global CONAN_APP
    CONAN_APP = shutil.which("conan")
    logger.info(f"Conan found: '{CONAN_APP}'")


def run_conan(args, **kwargs):
    if not "env" in kwargs:
        kwargs["env"] = CONAN_ENV
    else:
        kwargs["env"] |= CONAN_ENV
    kwargs["env"] |= os.environ
    kwargs["text"] = kwargs.get("text", True)
    args = [CONAN_APP] + args
    res = subprocess.run(args, shell=False, **kwargs)
    if res.returncode:
        logger.critical("Error while executing conan")
        print(res.stdout)
        print(res.stderr)
        exit(1)
    return res


def download(url, destdir):
    """Download file from url into destdir."""
    filename = destdir / "deps.zip"
    local_filename, _ = urlretrieve(url, filename=filename)
    downloaded = ZipFile(local_filename)
    downloaded.extractall(destdir)


def install(filename, destdir):
    """Install file from local directory into destdir."""
    logger.info(f"Importing {filename}")
    zipfile = ZipFile(str(filename))
    zipfile.extractall(destdir)


if __name__ == "__main__":

    # Set-up logger
    logger.setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO)

    # Get settings
    with open("luxcore.ini") as f:
        settings = configparser.ConfigParser(allow_no_value=True)
        settings.read_file(f)

    # Process
    with tempfile.TemporaryDirectory() as tmpdir:

        tmpdir = Path(tmpdir)

        CONAN_HOME = tmpdir / ".conan2"

        CONAN_ENV = {
            "CONAN_HOME": str(CONAN_HOME),
            "GCC_VERSION": str(settings["Build"]["gcc"]),
            "CXX_VERSION": str(settings["Build"]["cxx"]),
            "BUILD_TYPE": "Release",  # TODO Command line parameter
        }

        # Initialize
        ensure_conan_app()
        url = build_url(
            settings["Dependencies"]["user"],
            settings["Dependencies"]["release"]
        )
        gcc_version = settings["Build"]["gcc"]
        cxx_version = settings["Build"]["cxx"]

        # Download and unzip
        logger.info(f"Downloading dependencies (url='{url}')")
        download(url, tmpdir)

        # Clean
        logger.info("Cleaning local cache")
        res = run_conan(["remove", "-c", "*"], capture_output=True)
        for line in res.stderr.splitlines():
            logger.info(line)

        # Install
        logger.info("Installing")
        archive = tmpdir / "conan-cache-save.tgz"
        res = run_conan(
            ["cache", "restore", archive],
            capture_output=True,
        )
        for line in res.stderr.splitlines():
            logger.info(line)

        # Check
        logger.info("Checking integrity")
        res = run_conan(["cache", "check-integrity", "*"], capture_output=True)
        logger.info("Integrity check: OK")

        # Installing profiles
        logger.info("Installing profiles")
        run_conan(
            ["config", "install-pkg", "luxcoreconf/2.10.0@luxcore/luxcore"]
        )  # TODO version as a param

        # Generate & deploy
        logger.info("Generating")
        run_conan(
            [
                "install",
                "--requires=luxcoredeps/2.10.0@luxcore/luxcore",  # TODO version as a param
                "--build=missing",
                f"--profile:all={get_profile_name()}",
                "--deployer=full_deploy",
                f"--deployer-folder={BUILD_DIR}",
                "--generator=CMakeToolchain",
                "--generator=CMakeDeps",
                f"--output-folder={CMAKE_DIR}",
            ]
        )
