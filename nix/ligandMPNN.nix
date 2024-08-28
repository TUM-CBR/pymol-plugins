{
  packages,
  nixpkgs ? packages.nixpkgs,
  fetchzip ? nixpkgs.fetchzip,
  fetchurl ? nixpkgs.fetchurl,
  ligandMPNN ? fetchurl {
    url = "https://github.com/TUM-CBR/LigandMPNN/archive/refs/heads/main.zip";
    sha256 = "sha256-a1ixZutcaYVsueBajGKtEJlyKLnkdyVEL6DIYV5MuV4=";
    # sha256 = "sha256-joXI/JEgLcQb+LWzfUy9CMzUGCAzeYRmIDnWkSpvDrw=";
  }
}:
let
  inherit (nixpkgs) python311 python311Packages;
  model-urls = [
    {
      name = "ligandmpnn_v_32_005_25.pt";
      sha256 = "sha256-7genr7U7zpigsdM5lrvjtGsIMd9c4GhNPL1V7tGqkmM=";
    }
    {
      name = "ligandmpnn_v_32_010_25.pt";
      sha256 = "sha256-FhzSZAYf2paAy7lAJVUirkLylmxVLQRdh5E9lFKoCXA=";
    }
  ];
  fetchmpnn-model = { name, sha256 }: {
    name = name;
    model = fetchurl {
      url = "https://files.ipd.uw.edu/pub/ligandmpnn/${name}";
      sha256 = sha256;
    };
  };
  models = map fetchmpnn-model model-urls;
  biopython-179 = python311Packages.buildPythonPackage rec {
    pname = "biopython";
    version = "1.79";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/b6/a6/6e73bfa0a297b9884a029d224a8841b44b8dec35f17ceb1ef279d2a7dc76/biopython-1.79-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl";
      sha256 = "sha256-U1ynUGB4amguZXKr3EJCD6ilSvOIKX2nxW0VHHzGPuw=";
    };
    doCheck = false;
  };
  numpy-123 = python311Packages.buildPythonPackage rec {
    pname = "numpy";
    version = "1.23.5";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/e8/ad/b935c7421657a032fd2a5332eed098f3b9993a155afceb1daa280ff6611f/numpy-1.23.5-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl";
      sha256 = "sha256-WPVF79EQjmR2BKG1qoCVkczSVA9GiogL7blyR+cts4c=";
    };
  };
  pyparsing-311 = python311Packages.buildPythonPackage rec {
    pname = "pyparsing";
    version = "3.1.1";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/39/92/8486ede85fcc088f1b3dba4ce92dd29d126fd96b0008ea213167940a2475/pyparsing-3.1.1-py3-none-any.whl";
      sha256 = "sha256-MsfAtxFJPHL/GKmB0k8oqvnB+37V6WZ8noTj22I72/s=";
    };
  };
      
  prody = python311Packages.buildPythonPackage rec {
    pname = "ProDy";
    version = "2.4.1";
    src = python311Packages.fetchPypi {
      inherit pname version;
      sha256 = "sha256-nN+E2uNJs0Q3k7NVhjwhMWjMg32Wh0ziLoHu5OhSGcY=";
    };
    doCheck = false;
    propagatedBuildInputs = (with python311Packages; [ numpy-123 biopython-179 pyparsing-311 pip setuptools wheel ]);
  };
  lMPNNPython = python311.withPackages (ps: with ps; [ numpy-123 biopython-179 pyparsing-311 prody ]);
  mpnnModelsLinks =
    builtins.concatStringsSep "\n"
    (map ({name, model}: "ln -s ${model} $out/usr/share/ligandmpnn/model_params/${name}") models);
in
  nixpkgs.stdenv.mkDerivation rec {
    name = "ligandMPNN";
    src = ligandMPNN;
    buildInputs = [ lMPNNPython ];
    nativeBuildInputs = with nixpkgs; [ makeWrapper unzip ];
    mpnnModels = mpnnModelsLinks;
    configurePhase = "";
    buildPhase = "";
    installPhase = ''
      mkdir -p $out/bin
      mkdir -p $out/usr/share/ligandmpnn/model_params
      cp -r * $out/usr/share/ligandmpnn
      makeWrapper ${lMPNNPython}/bin/python $out/bin/ligandmpnn --add-flags $out/usr/share/ligandmpnn/run.py
      ${mpnnModelsLinks}
    '';
