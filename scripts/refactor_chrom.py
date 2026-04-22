import os
import re

def process_file(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    changed = False
    
    # 1. Replace struct fields and function signatures
    # Look for chrom: String or pub chrom: String
    new_content = re.sub(r'(\bchrom\s*:\s*)String\b', r'\1Arc<str>', content)
    if new_content != content:
        changed = True
        content = new_content

    if not changed:
        return False

    # 2. Add 'use std::sync::Arc;' if missing and Arc is now used
    if 'Arc<' in content and 'std::sync::Arc' not in content and 'use std::sync::Arc;' not in content:
        # Try to find a good place to add the import
        # Look for the first 'use ' line
        use_match = re.search(r'^use\s+', content, re.MULTILINE)
        if use_match:
            content = content[:use_match.start()] + "use std::sync::Arc;\n" + content[use_match.start():]
        else:
            # No 'use' lines, add after module doc comments if any
            doc_match = re.search(r'^(!|//)', content)
            if doc_match:
                # Find end of doc comments
                lines = content.split('\n')
                idx = 0
                for i, line in enumerate(lines):
                    if not line.startswith('//') and not line.startswith('!'):
                        idx = i
                        break
                lines.insert(idx, "use std::sync::Arc;")
                content = '\n'.join(lines)
            else:
                content = "use std::sync::Arc;\n" + content

    # 3. Simple fixes for String::new() or "".to_string() assigned to chrom
    # This is risky but can help
    # chrom: String::new() -> chrom: Arc::from("")
    content = re.sub(r'chrom\s*:\s*String::new\(\)', 'chrom: Arc::from("")', content)
    content = re.sub(r'chrom\s*:\s*"([^"]*)"\.to_string\(\)', r'chrom: Arc::from("\1")', content)
    
    # Fix common pattern: chrom: s.to_string() where s is &str
    # This is harder to do safely with regex

    with open(file_path, 'w') as f:
        f.write(content)
    
    return True

src_dir = 'src'
files_changed = 0
for root, dirs, files in os.walk(src_dir):
    for file in files:
        if file.endswith('.rs'):
            if process_file(os.path.join(root, file)):
                files_changed += 1
                print(f"Updated {os.path.join(root, file)}")

print(f"Total files updated: {files_changed}")
