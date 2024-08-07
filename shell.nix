{
  nixpkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.05.tar.gz") { },
  python3 ? nixpkgs.python3,
  fetchFromGitHub ? nixpkgs.fetchFromGitHub,
  nixGL ? import (fetchGit { url = "https://github.com/nix-community/nixGL.git"; rev = "310f8e49a149e4c9ea52f1adf70cdc768ec53f8a"; }) { },

}:
let
  pymol-src = python3.pkgs.buildPythonPackage {
    pname = "pymol-opensource-dev";
    version = "3.0.0";
    doCheck = false;
    src = fetchFromGitHub {
      owner = "netogallo";
      repo = "pymol-open-source";
      rev = "8cf28732dea7154db279c1d16704b7dfd65079f1";
      hash = "sha256-E6vd6nSoH+aTSDdY34aCDuhqKJSNKYjQaE40zJycXlg=";
    };
    buildInputs = with nixpkgs; [ python3.pkgs.numpy python3.pkgs.pyqt5 glew glm libpng libxml2 freetype msgpack netcdf ];
    env.NIX_CFLAGS_COMPILE = "-I ${nixpkgs.libxml2.dev}/include/libxml2";
    build-system = [
      python3.pkgs.setuptools
    ];
  };
  python-pymol-dev = python3.withPackages (
    py-pkgs: with py-pkgs; [
      pymol-src
      ipython
      pyqt5
      pyqt5.pyqt5-sip
      requests
      numpy
      biopython
      debugpy
      setuptools
      wheel
    ]
  );
in
nixpkgs.mkShell {
  name = "cbr-tools";
  packages = [
    python-pymol-dev
    nixGL
    nixpkgs.qt5.full
  ];

  shellHook = ''
    export PYMOL_STARTUP_PATH_POLY=$PWD
    export QT_QPA_PLATFORM=xcb
    unset SUDO_COMMAND
  '';
}

