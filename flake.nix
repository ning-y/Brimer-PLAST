{
  description = "Brimer-PLAST: local primer design with primer3-py + tnBLAST";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        # Custom tnBLAST derivation (minimal build, no MPI, no NCBI toolkit)
        tntblast = pkgs.stdenv.mkDerivation rec {
          pname = "tntblast";
          version = "2.77";

          src = pkgs.fetchFromGitHub {
            owner = "jgans";
            repo = "thermonucleotideBLAST";
            rev = "v${version}";
            hash = "sha256-PB3D2J5qXpnwvrrPjPkaoqQNLR2I6/915TVHVJzM+XM=";
          };

          nativeBuildInputs = [ pkgs.gnumake pkgs.gcc14 ];
          buildInputs = [ pkgs.zlib ];

          env.CC = "g++";
          env.CXX = "g++";

          patchPhase = ''
            echo "=== Patching Makefile for minimal tnBLAST build ==="
            # Comment out BLAST_DIR reference — Make's ifdef is true even for empty values
            sed -i Makefile \
              -e 's|^BLAST_DIR =.*|# BLAST_DIR commented out (minimal build)|' \
              -e 's|-DUSE_MPI||g' \
              -e 's|-DUSE_BLAST_DB||g' \
              -e 's|^CC = mpic++|CC = g++|'
            grep -n 'BLAST_DIR\|USE_MPI\|USE_BLAST\|CC =' Makefile | head -10
          '';

          buildPhase = ''
            make
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp tntblast $out/bin/
          '';

          meta = with pkgs.lib; {
            description = "DNA/RNA sequence database search with PCR primer and probe queries";
            homepage = "https://github.com/jgans/thermonucleotideBLAST";
            license = licenses.gpl2Only;
            platforms = platforms.linux;
          };
        };

        pythonEnv = pkgs.python312.withPackages (ps: with ps; [
          primer3
          typer
          pytest
          pyfaidx
          reportlab
        ]);

      in {
        packages = {
          inherit tntblast;
          brimer-plast = pythonEnv;
          default = pythonEnv;
        };

        devShells.default = pkgs.mkShell {
          packages = [
            pythonEnv
            tntblast
            pkgs.ruff
            pkgs.pyright
            pkgs.python312Packages.pip
          ];

          shellHook = ''
            # Create a venv if it doesn't exist, linking to the Nix Python
            if [ ! -d .venv ]; then
              ${pkgs.python312}/bin/python -m venv .venv --system-site-packages
            fi
            source .venv/bin/activate

            # Install the project in editable mode with dev extras
            if ! pip install -e ".[dev]" --quiet; then
              echo "Error: pip install failed. Check network connectivity or run 'pip install -e .' manually." >&2
              exit 1
            fi

            echo "Brimer-PLAST dev shell ready."
            echo "  Python: $(which python) ($(python --version))"
            echo "  tnBLAST: $(which tntblast || echo 'not built yet')"
          '';
        };
      });
}
