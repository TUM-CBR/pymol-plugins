{
  nixpkgs
}:
import (fetchGit {
  url = "https://github.com/TUM-CBR/cbr-tools-extra.git";
  rev = "6679418f50bb956119099465b43357e52cb20b5a";
}) { inherit nixpkgs; }
