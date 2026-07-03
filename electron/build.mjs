#!/usr/bin/env node
/**
 * Electron installer build script.
 *
 * Queries the Python package for the version (single source of truth from
 * setuptools-scm / git tag) and injects it into the electron-builder
 * metadata via --extraMetadata.version.
 *
 * Usage: node build.mjs <builder-flags...>
 *   node build.mjs --win
 *   node build.mjs --mac --arm64
 *   node build.mjs --linux AppImage
 */
import { execSync } from 'child_process';

// Get version from Python package (single source of truth).
let version;
try {
  version = execSync(
    'python3 -c "from brimer_plast import __version__; print(__version__)"',
    { encoding: 'utf-8', timeout: 10000 },
  ).trim();
} catch {
  version = '0.0.0';
}

// Check for bogus fallback versions before conversion.
// setuptools-scm's fallback_version is "0.0.0"; if no real git tag
// exists guess-next-dev bumps it to 0.0.1.devN+gXXXX.  Both are
// meaningless — use the sentinel so the UI filters them out.
if (/^0\.0\.\d/.test(version)) {
  version = '0.0.0';
} else {
  // Convert PEP 440 version to npm semver for electron-builder.
  // The only difference is the pre-release separator:
  //   PEP 440:  ".dev"  (e.g. 0.1.2.dev2+g9fc1cb508)
  //   npm:      "-dev." (e.g. 0.1.2-dev.2+g9fc1cb508)
  // Clean tags without a pre-release (0.1.1) pass through unchanged.
  version = version.replace(/\.dev(\d+)/, '-dev.$1');
}

const flagStr = process.argv.slice(2).join(' ');
execSync(`npx electron-builder ${flagStr} -c.extraMetadata.version=${version}`, {
  stdio: 'inherit',
  shell: true,
});