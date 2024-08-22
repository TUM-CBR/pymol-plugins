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
      wheel
    ]
  );
in
nixpkgs.mkShell {
  name = "cbr-tools";
  packages = with nixpkgs; [
    python-pymol-dev
    nixGL.nixGLMesa
    qt5.full
    qtcreator
    pyright
    ligandMPNN
  ];

  shellHook = ''
    export PYMOL_STARTUP_PATH_POLY=$PWD
    export QT_QPA_PLATFORM=xcb
    export DISPLAY=:1
    unset SUDO_COMMAND
  '';
}

