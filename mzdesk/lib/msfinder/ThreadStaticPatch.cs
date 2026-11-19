// mzDesk thread fix for MS-FINDER v3.73 StructureFinder.dll (compiled and run by `python -m mzdesk.ms2.msfinder_setup`).
//
// Why: MS-FINDER scores the structure candidates of a formula with Parallel.ForEach. FragmentPeakMatcher.GetSpectralAssignmentResult
// stores the largest bond dissociation energy of the candidate's fragments in the STATIC field maxBondDissociationEnergy
// (src/MSDIAL4/StructureFinderStandard/Spectral assigner/SpectralAssigner.cs, lines 21 and 38 at tag MSFINDER-v3.73), and
// getMatchedFragmentInfo reads it back (line 341) for the bond-energy likelihood. All threads share the one field, so a candidate is
// often scored with another candidate's value and the Total score changes from run to run (up to 0.5 on identical input).
//
// What: the field gets [ThreadStatic]. The value is written and read by the same thread within one call, so each thread keeps its own
// copy and every candidate is scored as the single-threaded code would score it. Nothing else in the assembly changes.
//
// Usage: mono ThreadStaticPatch.exe <StructureFinder.dll in> <StructureFinder.dll out>
// Exit 0 = written and verified; 2 = the field is already [ThreadStatic]; 1 = error.
using System;
using System.IO;
using System.Linq;
using Mono.Cecil;

static class ThreadStaticPatch
{
    const string TypeName = "Riken.Metabolomics.StructureFinder.SpectralAssigner.FragmentPeakMatcher";
    const string FieldName = "maxBondDissociationEnergy";
    const string Attribute = "System.ThreadStaticAttribute";

    static FieldDefinition FindField(AssemblyDefinition asm)
    {
        var type = asm.MainModule.GetType(TypeName) ?? throw new Exception("type " + TypeName + " not found");
        var field = type.Fields.FirstOrDefault(f => f.Name == FieldName) ?? throw new Exception("field " + FieldName + " not found");
        if (!field.IsStatic || field.FieldType.FullName != "System.Double") throw new Exception("field " + FieldName + " is not a static double");
        return field;
    }

    static bool IsThreadStatic(FieldDefinition f) => f.CustomAttributes.Any(a => a.AttributeType.FullName == Attribute);

    static int Main(string[] args)
    {
        try {
            if (args.Length != 2) { Console.Error.WriteLine("usage: ThreadStaticPatch.exe <in.dll> <out.dll>"); return 1; }
            var resolver = new DefaultAssemblyResolver();
            resolver.AddSearchDirectory(Path.GetDirectoryName(Path.GetFullPath(args[0])));
            var asm = AssemblyDefinition.ReadAssembly(args[0], new ReaderParameters { AssemblyResolver = resolver });
            var field = FindField(asm);
            if (IsThreadStatic(field)) { Console.WriteLine("already [ThreadStatic]: " + field.FullName); return 2; }

            // ThreadStaticAttribute through the assembly's own core reference (netstandard), so no new assembly reference is added.
            var module = asm.MainModule;
            var core = module.AssemblyReferences.FirstOrDefault(r => r.Name == "netstandard")
                       ?? module.AssemblyReferences.FirstOrDefault(r => r.Name == "mscorlib")
                       ?? throw new Exception("no netstandard / mscorlib reference");
            var attrType = new TypeReference("System", "ThreadStaticAttribute", module, core);
            var ctor = new MethodReference(".ctor", module.TypeSystem.Void, attrType) { HasThis = true };
            field.CustomAttributes.Add(new CustomAttribute(module.ImportReference(ctor)));
            asm.Write(args[1]);

            var check = FindField(AssemblyDefinition.ReadAssembly(args[1], new ReaderParameters { AssemblyResolver = resolver }));
            if (!IsThreadStatic(check)) throw new Exception("the written assembly does not carry the attribute");
            Console.WriteLine("patched: " + check.FullName + " is [ThreadStatic]");
            return 0;
        }
        catch (Exception e) {
            Console.Error.WriteLine("ThreadStaticPatch failed: " + e.Message);
            return 1;
        }
    }
}
