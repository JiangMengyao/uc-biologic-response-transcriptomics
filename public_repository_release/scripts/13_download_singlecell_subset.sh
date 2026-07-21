#!/usr/bin/env bash
set -euo pipefail

project_root="${PROJECT_ROOT:-$(pwd)}"
manifest_path="${project_root}/config/data_manifest.csv"
target_dir="${project_root}/data/raw/GSE282122"
mkdir -p "${target_dir}"
lock_dir="${target_dir}/.download.lock"
if ! mkdir "${lock_dir}" 2>/dev/null; then
  echo "Another single-cell download appears to be active: ${lock_dir}" >&2
  echo "If no downloader is running, remove this stale empty directory and retry." >&2
  exit 1
fi
trap 'rmdir "${lock_dir}" 2>/dev/null || true' EXIT

if [[ ! -f "${manifest_path}" ]]; then
  echo "Missing data manifest: ${manifest_path}" >&2
  exit 1
fi

required_files=(
  "paired_sample_list.csv"
  "UMAP_combined_objects.txt.gz"
  "myeloid_final.h5ad"
  "epicolonic_final.h5ad"
  "fibperi_final.h5ad"
)

file_md5() {
  local path="$1"
  if command -v md5 >/dev/null 2>&1; then
    md5 -q "${path}"
  elif command -v md5sum >/dev/null 2>&1; then
    md5sum "${path}" | awk '{print $1}'
  else
    echo "Neither md5 nor md5sum is available" >&2
    exit 1
  fi
}

manifest_field() {
  local filename="$1"
  local column="$2"
  awk -F',' -v file="${filename}" -v col="${column}" '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    $2==file { print $h[col]; exit }
  ' "${manifest_path}"
}

for filename in "${required_files[@]}"; do
  url="$(manifest_field "${filename}" source_url)"
  expected_md5="$(manifest_field "${filename}" expected_md5)"
  destination="${target_dir}/${filename}"
  if [[ -z "${url}" || -z "${expected_md5}" ]]; then
    echo "Manifest is incomplete for ${filename}" >&2
    exit 1
  fi

  if [[ -s "${destination}" ]]; then
    observed_md5="$(file_md5 "${destination}")"
    if [[ "${observed_md5}" == "${expected_md5}" ]]; then
      echo "CHECKSUM PASS ${filename}"
      continue
    fi
    echo "Checksum mismatch for existing ${filename}; resuming/replacing download" >&2
  fi

  curl --fail --location --retry 8 --retry-all-errors --retry-delay 3 \
    --connect-timeout 30 --speed-limit 1024 --speed-time 120 --continue-at - \
    --output "${destination}" "${url}"
  observed_md5="$(file_md5 "${destination}")"
  if [[ "${observed_md5}" != "${expected_md5}" ]]; then
    echo "CHECKSUM FAILED ${filename}: observed=${observed_md5} expected=${expected_md5}" >&2
    exit 1
  fi
  echo "CHECKSUM PASS ${filename}"
done

echo "All GSE282122 single-cell inputs are present and checksum-verified."
