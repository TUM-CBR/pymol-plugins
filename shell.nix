{
  packages ? import ../nix-packages/default.nix {},
  nixpkgs ? packages.nixpkgs,
  python3 ? nixpkgs.python3,
  nixGL ? packages.nixGL,
  fetchFromGitHub ? nixpkgs.fetchFromGitHub,
  ligandMPNN ? import ./nix/ligandMPNN.nix { inherit packages; }
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
      jedi
      jedi-language-server
      pyflakes
      pyqt5
      pyqt5.pyqt5-sip
      python-lsp-server
      requests
      numpy
      biopython
      debugpy
      setuptools
      virtualenv
      wheel
    ]
  );
in
nixpkgs.mkShell {
  name = "cbr-tools";
  packages = with nixpkgs; [
    clustal-omega
    python-pymol-dev
    #nixGL.nixGLMesa
    qt5.full
    qtcreator
    pyright
    libz
    gcc
    # ligandMPNN
  ];

  shellHook = ''
    export LD_LIBRARY_PATH=${nixpkgs.libz}/lib:/nix/store/p44qan69linp3ii0xrviypsw2j4qdcp2-gcc-13.2.0-lib/lib/:$LD_LIBRARY_PATH
    export PYMOL_STARTUP_PATH_POLY=$PWD
    export QT_QPA_PLATFORM=xcb
    unset SUDO_COMMAND
    export PATH=$PWD/tools/LigandMPNN:$PATH
    export THEPYTHON=${python-pymol-dev}
  '';
}

