{
  nixpkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.05.tar.gz") { },
  cbr-tools-extra ? import (fetchGit { url = "https://github.com/TUM-CBR/cbr-tools-extra.git"; rev = "fef22dcbd8880b2329f00c15c50ccaeca481aa92";}) { inherit nixpkgs; }, 
  nixGL ? import (fetchGit { url = "https://github.com/nix-community/nixGL.git"; rev = "310f8e49a149e4c9ea52f1adf70cdc768ec53f8a"; }) { },
  pymol ? nixpkgs.pymol,
  lib ? nixpkgs.lib
}:
let
  python = lib.findFirst
    (p: p?pythonAtLeast && p.pythonAtLeast "3")
    (throw "Python was expected in the 'nativeBuildInputs' of 'pymol'.")
    pymol-ext.nativeBuildInputs;
  cbr-tools =
    nixpkgs.stdenv.mkDerivation {
      name = "cbr-tools";
      nativeBuildInputs = [
	nixpkgs.coreutils
      ];
      buildInputs = [
	cbr-tools-extra
      ];
      src = ./cbr;
      installPhase = ''
	plugin_folder="$out/${python.sitePackages}/pymol/pymol_path/data/startup/"
	mkdir -p $plugin_folder
	cp -r "$PWD" $plugin_folder/
      '';
    };

  # The pymol distribution from schrodinger bundles more python packages
  # than the pymol available in nixpkgs
  pymol-required-pkgs = python.pkgs.requiredPythonModules (with python.pkgs; [ pyqt5 pyqt5.pyqt5-sip requests numpy ]); 
  pymol-ext = pymol.overridePythonAttrs {
    buildInputs = pymol.buildInputs ++ (with python.pkgs; [ requests ]);
    doCheck = false;
    postInstall = ''
      wrapProgram $out/bin/pymol \
	--prefix PYTHONPATH : ${lib.makeSearchPathOutput "lib" python.sitePackages pymol-required-pkgs}
    '';
  };
  pymolpath-out = "${placeholder "out"}/${python.sitePackages}/pymol/pymol_path/";
  pymol-exe = "${pymol-ext}/bin/pymol";
  pymol-launcher = ''#!/bin/bash
  ${nixGL}/bin/nixGL ${placeholder "out"}/bin/pymol-unwrapped
  '';
in
nixpkgs.buildEnv {
  name = "pymol-cbr-tools";
  paths = [ pymol-ext nixGL cbr-tools python python.pkgs.requests ];
  nativeBuildInputs = [ nixpkgs.makeBinaryWrapper ];
  postBuild = ''
    if [ -L $out/bin ]; then
      unlink "$out/bin"
    fi

    mkdir -p $out/bin

    makeWrapper ${pymol-exe} "$out/bin/pymol-unwrapped" --set PYMOL_PATH "${pymolpath-out}" --set QT_QPA_PLATFORM "xcb" 
    pymol_wrapper=$out/bin/pymol
    '';
}

