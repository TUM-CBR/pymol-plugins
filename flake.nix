{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-24.05";
    flake-parts.url = "github:hercules-ci/flake-parts";
    cbr-tools-extra-flake = {
      url = "github:TUM-CBR/cbr-tools-extra?ref=fef22dcbd8880b2329f00c15c50ccaeca481aa92";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    nixGL-flake = {
      url = "github:nix-community/nixGL?ref=310f8e49a149e4c9ea52f1adf70cdc768ec53f8a";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = inputs@{ self, nixpkgs, flake-parts, cbr-tools-extra-flake, nixGL-flake, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
	"x86_64-linux"
	"x86_64-darwin"
      ];
      perSystem = { self', pkgs, system, ... }:
      let
	nixGL = nixGL-flake.outputs.packages.${system}.default;
	cbr-tools-extra = cbr-tools-extra-flake.outputs.packages.${system}.default;
	cbr-tools = import ./default.nix {
	  nixpkgs = pkgs;
	  inherit nixGL cbr-tools-extra;
	};
      in
      {
	packages = {
	  inherit cbr-tools;
	  default = cbr-tools;
	};
      };
    };
}
