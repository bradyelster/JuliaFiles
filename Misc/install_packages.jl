using Pkg

# Define the path to the text file
package_file = "juliareqs.txt"

# Read the file and split into lines (package names)
package_names = readlines(package_file)

# Loop over each package name and install it
for pkg in package_names
    Pkg.add(pkg)
end
