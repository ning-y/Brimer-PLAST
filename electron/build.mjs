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

// Strip everything after the major.minor.patch triplet so the result
// is always valid npm semver for electron-builder.  Pre-release suffixes
// (e.g. .dev3) and local suffixes (+gdeadbeef) are dropped.
version = version.replace(/^(\d+\.\d+\.\d+).*$/, '$1');

const flagStr = process.argv.slice(2).join(' ');
execSync(`npx electron-builder ${flagStr} -c.extraMetadata.version=${version}`, {
  stdio: 'inherit',
  shell: true,
});