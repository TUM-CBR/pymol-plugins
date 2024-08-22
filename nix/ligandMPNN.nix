{
  packages,
  nixpkgs ? packages.nixpkgs,
  fetchzip ? nixpkgs.fetchzip,
  ligandMPNN ? fetchzip {
    url = "https://github.com/TUM-CBR/LigandMPNN/archive/refs/heads/main.zip";
    sha256 = "sha256-joXI/JEgLcQb+LWzfUy9CMzUGCAzeYRmIDnWkSpvDrw=";
  }
}:
let
  inherit (nixpkgs) python311 python311Packages;
  biopython-179 = python311Packages.biopython.overridePythonAttrs {
    src = fetchzip {
      url = "https://github.com/biopython/biopython/archive/refs/tags/biopython-179.tar.gz";
      sha256 = "sha256-nBqvfVbuv2yVgVc451xRRwt2dycQZzWrWDP0GuwEkfc=";
    };
    doCheck = false;
  };
  prody = python311Packages.buildPythonPackage rec {
    pname = "ProDy";
    version = "2.4.1";
    src = python311Packages.fetchPypi {
      inherit pname version;
      sha256 = "sha256-nN+E2uNJs0Q3k7NVhjwhMWjMg32Wh0ziLoHu5OhSGcY=";
    };
    doCheck = false;
    propagatedBuildInputs = with python311Packages; [ biopython-179 pyparsing pip setuptools numpy scipy ];
  };
  lMPNNPython = python311.withPackages (ps: with ps; [ torch prody ]);
in
  nixpkgs.runCommand "ligandMPNN" { buildInputs = [nixpkgs.makeBinaryWrapper]; } ''
    mkdir -p $out/bin
    makeWrapper ${lMPNNPython}/bin/python $out/bin/ligandmpnn --add-flags ${ligandMPNN}/run.py
    ''
