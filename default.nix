{
  nixpkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.05.tar.gz") { },
  cbr-tools-extra ? import (fetchGit { url = "https://github.com/TUM-CBR/cbr-tools-extra.git"; rev = "fef22dcbd8880b2329f00c15c50ccaeca481aa92";}) { inherit nixpkgs; }, 
  nixGL ? import (fetchGit { url = "https://github.com/nix-community/nixGL.git"; rev = "310f8e49a149e4c9ea52f1adf70cdc768ec53f8a"; }) { },
  pymol ? nixpkgs.pymol,
  lib ? nixpkgs.lib,
  raspp-src ? (fetchGit { url = "https://github.com/TUM-CBR/SCHEMA-RASPP.git"; rev="bd3befee0e99b97c28200b385f242e2e3abf9ec6"; })
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
      src = ./cbr;
      inherit raspp-src;
      installPhase = ''
	plugin_folder="$out/${python.sitePackages}/pymol/pymol_path/data/startup/"
	mkdir -p $plugin_folder
	cp -r "$PWD" $plugin_folder/
	raspp_folder="$plugin_folder/cbr/schema/raspp"
	cp -r "${raspp-src}" $raspp_folder 
	ls $plugin_folder/cbr/schema/raspp
      '';
    };

  # The pymol distribution from schrodinger bundles more python packages
  # than the pymol available in nixpkgs
  pymol-extra-deps = with python.pkgs; [
    pyqt5
    pyqt5-sip
    requests
    numpy
    biopython
  ];
  pymol-required-pkgs = python.pkgs.requiredPythonModules pymol-extra-deps; 
  pymol-ext = pymol.overridePythonAttrs {
    buildInputs = pymol.buildInputs ++ pymol-extra-deps;
    doCheck = false;
    postInstall = ''
      wrapProgram $out/bin/pymol \
	--prefix PYTHONPATH : ${lib.makeSearchPathOutput "lib" python.sitePackages pymol-required-pkgs}
    '';
  };
  pymolpath-out = "${placeholder "out"}/${python.sitePackages}/pymol/pymol_path/";
  pymol-exe = "${pymol-ext}/bin/pymol";
  nixGL-exe = "${nixGL}/bin/nixGL";
in
nixpkgs.buildEnv {
  name = "pymol-cbr-tools";
  paths = [ pymol-ext nixGL cbr-tools cbr-tools-extra ];
  nativeBuildInputs = [ nixpkgs.makeBinaryWrapper ];
  postBuild = ''
    if [ -L $out/bin ]; then
      unlink "$out/bin"
    fi

    mkdir -p $out/bin
    makeWrapper ${nixGL-exe} "$out/bin/pymol" --add-flags ${pymol-exe} --set PYMOL_PATH "${pymolpath-out}" --set QT_QPA_PLATFORM "xcb" 
    '';
}

