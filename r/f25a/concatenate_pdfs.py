#!/usr/bin/env python3
"""
Script to concatenate PDFs based on fstem values from YAML files.
The sorting order is: Fig < ExtFig < SuppFig (RevFig ignored)
Within each prefix: sort by integer, then letter, then second integer.
"""

import glob
import re
import subprocess
from pathlib import Path
from typing import List, Tuple

import yaml


def parse_fstem(fstem: str) -> Tuple[str, int, str, int]:
    """
    Parse fstem according to pattern: "<prefix><first positive integer><optional single letter><optional second positive integer>"

    Returns: (prefix, first_int, letter, second_int)
    Where letter is empty string if not present, second_int is 0 if not present
    """
    # Pattern: (prefix)(number)(optional letter)(optional number)
    pattern = r"^(Fig|ExtFig|SuppFig|RevFig)(\d+)([A-Z]?)(\d*)$"
    match = re.match(pattern, fstem)

    if not match:
        raise ValueError(f"Invalid fstem format: {fstem}")

    prefix = match.group(1)
    first_int = int(match.group(2))
    letter = match.group(3) if match.group(3) else ""
    second_int = int(match.group(4)) if match.group(4) else 0

    return (prefix, first_int, letter, second_int)


def sort_key(fstem: str) -> Tuple:
    """
    Create sort key for fstem according to the rules:
    1. Fig < ExtFig < SuppFig (RevFig ignored)
    2. Sort by first integer (ascending)
    3. Sort by letter (A < B < C...)
    4. Sort by second integer (ascending)
    """
    prefix, first_int, letter, second_int = parse_fstem(fstem)

    # Prefix priority: Fig=0, ExtFig=1, SuppFig=2, RevFig=999 (to be filtered out)
    prefix_priority = {
        "Fig": 0,
        "ExtFig": 1,
        "SuppFig": 2,
        "RevFig": 999,  # Will be filtered out
    }

    return (prefix_priority[prefix], first_int, letter, second_int)


def extract_fstem_from_yaml(yaml_file: str) -> List[str]:
    """Extract all fstem values from a YAML file."""
    fstems = []
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

        if isinstance(data, list):
            for item in data:
                if isinstance(item, dict) and "fstem" in item:
                    fstems.append(item["fstem"])

    return fstems


def find_all_fstems() -> List[str]:
    """Find all fstem values from all YAML files in notebooks directory."""
    all_fstems = []

    # Find all *_params.yml files
    yaml_files = glob.glob("notebooks/**/*_params.yml", recursive=True)

    for yaml_file in yaml_files:
        fstems = extract_fstem_from_yaml(yaml_file)
        all_fstems.extend(fstems)
        print(f"Found {len(fstems)} fstems in {yaml_file}")

    return all_fstems


def find_pdf_files() -> List[str]:
    """Find all PDF files in _output directory."""
    pdf_files = glob.glob("_output/**/*.pdf", recursive=True)
    return pdf_files


def main():
    print("=== PDF Concatenation Script ===")

    # Step 1: Extract all fstem values
    print("Step 1: Extracting fstem values from YAML files...")
    all_fstems = find_all_fstems()
    print(f"Total fstems found: {len(all_fstems)}")

    # Step 2: Filter out RevFig entries and sort
    print("Step 2: Filtering and sorting fstems...")
    valid_fstems = [fstem for fstem in all_fstems if not fstem.startswith("RevFig")]
    print(f"Valid fstems (excluding RevFig): {len(valid_fstems)}")

    # Sort according to rules
    sorted_fstems = sorted(valid_fstems, key=sort_key)

    print("Sorted fstems:")
    for i, fstem in enumerate(sorted_fstems, 1):
        print(f"  {i:2d}. {fstem}")

    # Step 3: Find corresponding PDF files
    print("Step 3: Finding corresponding PDF files...")
    pdf_files = find_pdf_files()
    pdf_dict = {}

    # Create a mapping from filename to full path
    for pdf_path in pdf_files:
        filename = Path(pdf_path).stem  # Get filename without extension
        pdf_dict[filename] = pdf_path

    # Find PDFs for each fstem
    found_pdfs = []
    missing_pdfs = []

    for fstem in sorted_fstems:
        if fstem in pdf_dict:
            found_pdfs.append(pdf_dict[fstem])
            print(f"  ✓ Found: {fstem} -> {pdf_dict[fstem]}")
        else:
            missing_pdfs.append(fstem)
            print(f"  ✗ Missing: {fstem}")

    print(f"\nFound {len(found_pdfs)} PDFs, missing {len(missing_pdfs)} PDFs")

    if missing_pdfs:
        print("Missing PDFs:")
        for missing in missing_pdfs:
            print(f"  - {missing}")

    # Step 4: Concatenate PDFs
    if found_pdfs:
        print("Step 4: Concatenating PDFs...")
        output_file = "concatenated_supplements.pdf"

        # Use ghostscript to concatenate PDFs
        cmd = [
            "gs",
            "-dBATCH",
            "-dNOPAUSE",
            "-q",
            "-sDEVICE=pdfwrite",
            f"-sOutputFile={output_file}",
        ] + found_pdfs

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"✓ Successfully created {output_file}")
            print(f"  Combined {len(found_pdfs)} PDFs in the correct order")
        except subprocess.CalledProcessError as e:
            print(f"✗ Error concatenating PDFs: {e}")
            print(f"  stderr: {e.stderr}")

            # Try alternative with pypdf if ghostscript fails
            try:
                print("Trying alternative method with pypdf...")
                import pypdf

                writer = pypdf.PdfWriter()

                for pdf_path in found_pdfs:
                    reader = pypdf.PdfReader(pdf_path)
                    for page in reader.pages:
                        writer.add_page(page)

                with open(output_file, "wb") as output:
                    writer.write(output)

                print(f"✓ Successfully created {output_file} using pypdf")
            except ImportError:
                print("pypdf not available. Please install ghostscript or pypdf2")
            except Exception as e:
                print(f"✗ Error with pypdf: {e}")
    else:
        print("No PDFs found to concatenate!")


if __name__ == "__main__":
    main()
